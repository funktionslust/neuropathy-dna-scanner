//! `MemoryBamSource` - an `AlignmentSource` over an in-memory BAM byte
//! buffer. Available on every target, including `wasm32`.
//!
//! Strategy: on every `query_region` call, construct a fresh
//! `noodles_bam::io::Reader` over an `std::io::Cursor<&[u8]>` and iterate
//! all records, filtering by chromosome + coordinate in-process. This is
//! deliberately naive - for each query we re-read the whole BAM.
//!
//! Why naive: the BAM index (BAI) is an on-disk structure that `noodles-bam`
//! currently only loads via a path-based builder. Teaching it to load from
//! in-memory BAI bytes is a follow-up. For a
//! 119 MB BAM slice at three regions the naive cost is fine - ~3 full scans,
//! which the POC measures.
//!
//! This source is used by:
//! - the WASM entry point (`wasm/src/lib.rs::analyze_buffer`),
//! - any core unit test that wants to exercise the real noodles decode path
//!   without touching the filesystem.

use std::io::Cursor;

use crate::reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region};
use crate::sources::{convert_trait_record, sam_header_to_alignment_header};
use crate::{Error, Result};

/// In-memory BAM alignment source.
///
/// Holds the full BAM bytes and (optionally) a cached parsed header. Each
/// `query_region` call walks all records in the buffer and filters by the
/// caller's region. Header parsing is cached after the first call.
pub struct MemoryBamSource {
    bytes: Vec<u8>,
    cached_header: Option<AlignmentHeader>,
    cached_noodles_header: Option<noodles_sam::Header>,
}

impl MemoryBamSource {
    /// Construct a new source from owned BAM bytes. The bytes are moved
    /// into the source and never mutated.
    pub fn new(bytes: Vec<u8>) -> Self {
        Self {
            bytes,
            cached_header: None,
            cached_noodles_header: None,
        }
    }

    fn parse_header(&mut self) -> Result<()> {
        if self.cached_header.is_some() {
            return Ok(());
        }
        let mut reader = noodles_bam::io::Reader::new(Cursor::new(&self.bytes));
        let h = reader
            .read_header()
            .map_err(|e| Error::Io(std::io::Error::other(format!("read BAM header: {e}"))))?;
        let converted = sam_header_to_alignment_header(&h);
        self.cached_header = Some(converted);
        self.cached_noodles_header = Some(h);
        Ok(())
    }

    fn noodles_header(&self) -> &noodles_sam::Header {
        self.cached_noodles_header
            .as_ref()
            .expect("parse_header called first")
    }
}

impl AlignmentSource for MemoryBamSource {
    fn read_header(&mut self) -> Result<AlignmentHeader> {
        self.parse_header()?;
        Ok(self.cached_header.as_ref().unwrap().clone())
    }

    fn query_region<'a>(
        &'a mut self,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<AlignmentRecord>> + 'a>> {
        self.parse_header()?;
        let target_ref_id: Option<usize> = self
            .noodles_header()
            .reference_sequences()
            .get_index_of(region.chrom.as_bytes());
        let Some(target_ref_id) = target_ref_id else {
            return Ok(Box::new(std::iter::empty()));
        };

        let region_start = region.start;
        let region_end = region.end;

        // Rebuild the reader from the start to iterate all records. We
        // accumulate matches into an owned Vec so the returned iterator
        // does not need to borrow the reader or the header across calls.
        let mut reader = noodles_bam::io::Reader::new(Cursor::new(&self.bytes));
        reader
            .read_header()
            .map_err(|e| Error::Io(std::io::Error::other(format!("skip BAM header: {e}"))))?;

        let mut collected: Vec<Result<AlignmentRecord>> = Vec::new();
        for rec_result in reader.records() {
            match rec_result {
                Ok(record) => {
                    let Some(Ok(ref_id)) = record.reference_sequence_id() else {
                        continue;
                    };
                    if ref_id != target_ref_id {
                        continue;
                    }
                    let Some(arec) = convert_trait_record(&record) else {
                        continue;
                    };
                    // Region overlap on 1-based inclusive coordinates.
                    if arec.end < region_start || arec.start > region_end {
                        continue;
                    }
                    collected.push(Ok(arec));
                }
                Err(e) => collected.push(Err(Error::Io(std::io::Error::other(format!(
                    "BAM record decode: {e}"
                ))))),
            }
        }
        Ok(Box::new(collected.into_iter()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Read as _;

    /// Reuse the minimal BAM fixture from the file source's test module by
    /// writing one to a tempfile, then slurping the bytes into memory.
    fn build_minimal_bam_bytes() -> Vec<u8> {
        use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().to_path_buf();
        {
            const MD5_HEX: &str = "abcdef0123456789abcdef0123456789";
            let mut ref_seq: Map<ReferenceSequence> =
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(83_257_441usize).unwrap());
            let tag =
                noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
            ref_seq.other_fields_mut().insert(tag, MD5_HEX.into());
            let header = noodles_sam::Header::builder()
                .add_reference_sequence("chr17", ref_seq)
                .build();
            let file = File::create(&path).unwrap();
            let mut writer = noodles_bam::io::Writer::new(file);
            writer.write_header(&header).unwrap();
            writer.try_finish().unwrap();
        }
        let mut bytes = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut bytes).unwrap();
        bytes
    }

    #[test]
    fn memory_bam_source_reads_header() {
        let bytes = build_minimal_bam_bytes();
        let mut src = MemoryBamSource::new(bytes);
        let header = src.read_header().unwrap();
        assert_eq!(header.sequences.len(), 1);
        assert_eq!(header.sequences[0].name, "chr17");
        assert_eq!(header.sequences[0].length, 83_257_441);
    }

    #[test]
    fn memory_bam_source_empty_region_yields_no_records() {
        let bytes = build_minimal_bam_bytes();
        let mut src = MemoryBamSource::new(bytes);
        let region = Region::new("chr17", 1, 1_000);
        let iter = src.query_region(&region).unwrap();
        assert_eq!(iter.count(), 0);
    }

    #[test]
    fn memory_bam_source_unknown_chrom_yields_no_records() {
        let bytes = build_minimal_bam_bytes();
        let mut src = MemoryBamSource::new(bytes);
        let region = Region::new("chrZZ", 1, 1_000);
        let iter = src.query_region(&region).unwrap();
        assert_eq!(iter.count(), 0);
    }
}
