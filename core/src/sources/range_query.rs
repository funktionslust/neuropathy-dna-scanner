//! BAI-driven byte range computation for sparse BAM fetching.
//!
//! Given a BAI index, a BAM header (for chromosome name → id lookup),
//! and a list of query regions, returns the file byte ranges the caller
//! must fetch to cover all overlapping BGZF blocks.
//!
//! Used by the WASM path: JS calls this to learn which
//! byte ranges to `File.slice()` from a large BAM, then passes those
//! ranges to `CursorBamSource` via an `OffsetCursor`.

use std::io::Cursor;

use noodles_csi::binning_index::BinningIndex;

use crate::reader::Region;
use crate::{Error, Result};

/// Maximum size of a compressed BGZF block (64 KB). Used as a safety
/// margin when computing the end of a chunk's byte range.
const MAX_BGZF_BLOCK_SIZE: u64 = 65536;

/// A file byte range: `(offset, length)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ByteRange {
    pub offset: u64,
    pub length: u64,
}

/// Compute the file byte ranges needed to decode all BAM records overlapping
/// the given regions.
///
/// Returns a sorted, merged list of `ByteRange` values. The first range
/// always covers the BAM header (offset 0).
///
/// # Arguments
/// - `bai_bytes`: the full BAI index file bytes.
/// - `header_bytes`: the BAM file header (first N bytes - typically <=64 KB).
///   Only used to map chromosome names to reference sequence IDs. The length
///   of this slice is also used as the header range's size.
/// - `regions`: the genomic regions to query.
pub fn required_ranges_bam(
    bai_bytes: &[u8],
    header_bytes: &[u8],
    regions: &[Region],
) -> Result<Vec<ByteRange>> {
    // Parse the BAI.
    let index: noodles_bam::bai::Index = noodles_bam::bai::io::Reader::new(Cursor::new(bai_bytes))
        .read_index()
        .map_err(|e| Error::Internal(format!("failed to parse BAI: {e}")))?;

    // Parse the BAM header to get the reference sequence dictionary.
    let mut bam_reader = noodles_bam::io::Reader::new(Cursor::new(header_bytes));
    let sam_header = bam_reader
        .read_header()
        .map_err(|e| Error::Internal(format!("failed to parse BAM header: {e}")))?;

    let mut raw_ranges: Vec<(u64, u64)> = Vec::new();

    // Always include the header.
    raw_ranges.push((0, header_bytes.len() as u64));

    for region in regions {
        let ref_id = match sam_header
            .reference_sequences()
            .get_index_of(region.chrom.as_bytes())
        {
            Some(id) => id,
            None => continue, // skip chromosomes not in the file
        };

        let start = noodles_core::Position::try_from(region.start as usize)
            .map_err(|_| Error::Internal(format!("invalid start: {}", region.start)))?;
        let end = noodles_core::Position::try_from(region.end as usize)
            .map_err(|_| Error::Internal(format!("invalid end: {}", region.end)))?;
        let interval = noodles_core::region::Interval::from(start..=end);

        let chunks = index
            .query(ref_id, interval)
            .map_err(|e| Error::Internal(format!("BAI query failed: {e}")))?;

        for chunk in &chunks {
            let start_offset = chunk.start().compressed();
            let end_offset = chunk.end().compressed();
            // Conservative: include the full BGZF block at the end.
            let range_end = end_offset + MAX_BGZF_BLOCK_SIZE;
            raw_ranges.push((start_offset, range_end - start_offset));
        }
    }

    Ok(merge_ranges(&mut raw_ranges))
}

/// Compute file byte ranges for CRAM containers covering the given regions.
///
/// CRAI records are a flat list (simpler than BAI's binning index). Each
/// record has `offset()` (container start), `landmark()` (slice offset in
/// container), `slice_length()` (compressed slice size), plus alignment
/// start/span for region overlap checking.
/// Default position margin for SV analysis (large regions, few queries).
pub const CRAM_MARGIN_DEFAULT: u32 = 5_000_000;
/// Tight position margin for SNV pileup (many small single-position queries).
pub const CRAM_MARGIN_TIGHT: u32 = 50_000;

pub fn required_ranges_cram(
    crai_bytes: &[u8],
    header_bytes: &[u8],
    regions: &[Region],
) -> Result<Vec<ByteRange>> {
    required_ranges_cram_with_margin(crai_bytes, header_bytes, regions, CRAM_MARGIN_DEFAULT)
}

pub fn required_ranges_cram_with_margin(
    crai_bytes: &[u8],
    header_bytes: &[u8],
    regions: &[Region],
    margin: u32,
) -> Result<Vec<ByteRange>> {
    let index: Vec<noodles_cram::crai::Record> =
        noodles_cram::crai::io::Reader::new(Cursor::new(crai_bytes))
            .read_index()
            .map_err(|e| Error::Internal(format!("failed to parse CRAI: {e}")))?;

    // Parse the CRAM header to get the reference sequence dictionary.
    let mut cram_reader =
        noodles_cram::io::reader::Builder::default().build_from_reader(Cursor::new(header_bytes));
    let sam_header = cram_reader
        .read_header()
        .map_err(|e| Error::Internal(format!("failed to parse CRAM header: {e}")))?;

    let mut raw_ranges: Vec<(u64, u64)> = Vec::new();

    // Include the full CRAM header. The header container may be larger
    // than the header_bytes slice we received - use the minimum data
    // container offset from the CRAI as the header end.
    let first_data_offset = index
        .iter()
        .filter(|r| r.reference_sequence_id().is_some())
        .map(|r| r.offset())
        .min()
        .unwrap_or(0);
    let header_end = if first_data_offset > 0 {
        first_data_offset
    } else {
        header_bytes.len() as u64
    };
    raw_ranges.push((0, header_end));

    // Per-region strategy: for each query region, find the CRAI
    // containers that overlap [region.start - margin, region.end + margin].
    // This produces tight byte ranges even when regions are scattered
    // across a chromosome. The noodles skip-before patch handles
    // decompression context without needing data from chromosome start.
    for region in regions {
        let ref_id = match sam_header
            .reference_sequences()
            .get_index_of(region.chrom.as_bytes())
        {
            Some(id) => id,
            None => continue,
        };

        let min_pos = region.start.saturating_sub(margin);
        let max_pos = region.end + margin;
        let mut min_offset: Option<u64> = None;
        let mut max_end: Option<u64> = None;

        for record in &index {
            let rec_ref_id = match record.reference_sequence_id() {
                Some(id) => id,
                None => continue,
            };
            if rec_ref_id != ref_id {
                continue;
            }

            let rec_start = record
                .alignment_start()
                .map(|p| p.get() as u32)
                .unwrap_or(0);
            let rec_span = record.alignment_span() as u32;
            let rec_end = rec_start.saturating_add(rec_span);

            if rec_end < min_pos || rec_start > max_pos {
                continue;
            }

            let container_start = record.offset();
            let container_end = container_start + record.landmark() + record.slice_length();

            min_offset = Some(min_offset.map_or(container_start, |m: u64| m.min(container_start)));
            max_end = Some(max_end.map_or(container_end, |m: u64| m.max(container_end)));
        }

        if let (Some(start), Some(end)) = (min_offset, max_end) {
            raw_ranges.push((start, end - start + 65536));
        }
    }

    Ok(merge_ranges(&mut raw_ranges))
}

/// Sort and merge overlapping/adjacent ranges into a minimal list.
fn merge_ranges(ranges: &mut [(u64, u64)]) -> Vec<ByteRange> {
    ranges.sort_by_key(|&(offset, _)| offset);

    let mut merged: Vec<ByteRange> = Vec::new();
    for &(offset, length) in ranges.iter() {
        let end = offset + length;
        if let Some(last) = merged.last_mut() {
            let last_end = last.offset + last.length;
            if offset <= last_end {
                // Overlapping or adjacent - extend.
                last.length = (end.max(last_end)) - last.offset;
                continue;
            }
        }
        merged.push(ByteRange { offset, length });
    }
    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Read as _;
    use std::path::PathBuf;

    fn fixture_path(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data")
            .join(name)
    }

    fn load(name: &str) -> Vec<u8> {
        let mut bytes = Vec::new();
        File::open(fixture_path(name))
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();
        bytes
    }

    #[test]
    fn required_ranges_bam_returns_non_empty_for_fixture() {
        let bai_bytes = load("cn3_30x.bam.bai");
        // Read just the first 64 KB as the "header bytes".
        let full_bam = load("cn3_30x.bam");
        let header_len = full_bam.len().min(65536);
        let header_bytes = &full_bam[..header_len];

        let regions = vec![
            Region::new("chr17", 15_229_777, 15_265_079), // PMP22
            Region::new("chr2", 50_000_000, 50_100_000),  // control
        ];

        let ranges = required_ranges_bam(&bai_bytes, header_bytes, &regions).unwrap();

        assert!(
            !ranges.is_empty(),
            "should return at least the header range"
        );
        // First range should start at 0 (header).
        assert_eq!(ranges[0].offset, 0);
        // Ranges should be sorted.
        for w in ranges.windows(2) {
            assert!(
                w[0].offset + w[0].length <= w[1].offset,
                "ranges should be sorted and non-overlapping"
            );
        }
        // Total bytes should be less than the full BAM.
        let total: u64 = ranges.iter().map(|r| r.length).sum();
        assert!(
            total < full_bam.len() as u64,
            "sparse ranges ({total} bytes) should be smaller than full BAM ({} bytes)",
            full_bam.len()
        );
    }

    #[test]
    fn merge_ranges_merges_overlapping() {
        let mut raw = vec![(0, 100), (50, 100), (200, 50)];
        let merged = merge_ranges(&mut raw);
        assert_eq!(merged.len(), 2);
        assert_eq!(
            merged[0],
            ByteRange {
                offset: 0,
                length: 150
            }
        );
        assert_eq!(
            merged[1],
            ByteRange {
                offset: 200,
                length: 50
            }
        );
    }

    #[test]
    fn merge_ranges_handles_adjacent() {
        let mut raw = vec![(0, 100), (100, 50)];
        let merged = merge_ranges(&mut raw);
        assert_eq!(merged.len(), 1);
        assert_eq!(
            merged[0],
            ByteRange {
                offset: 0,
                length: 150
            }
        );
    }
}
