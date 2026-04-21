//! Concrete [`crate::AlignmentSource`] implementations.
//!
//! - [`file`] - native `FileAlignmentSource` over `std::fs::File`. Gated off
//!   on `wasm32` targets because it needs the filesystem.
//! - [`memory`] - `MemoryBamSource` over an in-memory BAM byte buffer.
//!   Available on every target, including `wasm32`.
//! - [`offset_cursor`] - `OffsetCursor`, a sparse `Read + Seek` over
//!   prefetched byte ranges at known file offsets. Used by the WASM path.
//! - [`cursor_source`] - `CursorBamSource`, an `AlignmentSource` wrapping
//!   `OffsetCursor` in a noodles `IndexedReader`.

pub mod cursor_cram_source;
pub mod cursor_source;
#[cfg(not(target_arch = "wasm32"))]
pub mod file;
pub mod memory;
pub mod offset_cursor;
pub mod range_query;

use crate::reader::{AlignmentHeader, AlignmentRecord, SequenceInfo};

/// Convert a noodles SAM header into our stripped-down [`AlignmentHeader`]
/// shape. Shared by `FileAlignmentSource` and `MemoryBamSource` - both
/// decode BAM headers identically, and neither touches the filesystem in
/// this step.
pub fn sam_header_to_alignment_header(header: &noodles_sam::Header) -> AlignmentHeader {
    use noodles_sam::header::record::value::map::reference_sequence::tag;

    let sequences: Vec<SequenceInfo> = header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| {
            let md5 = seq
                .other_fields()
                .get(&tag::MD5_CHECKSUM)
                .and_then(|v| hex_to_md5_bytes(v.as_ref()));
            SequenceInfo {
                name: name.to_string(),
                length: usize::from(seq.length()) as u64,
                md5,
            }
        })
        .collect();
    AlignmentHeader { sequences }
}

pub(crate) fn hex_to_md5_bytes(hex: &[u8]) -> Option<[u8; 16]> {
    if hex.len() != 32 {
        return None;
    }
    let mut out = [0u8; 16];
    for (i, pair) in hex.chunks_exact(2).enumerate() {
        let hi = hex_nibble(pair[0])?;
        let lo = hex_nibble(pair[1])?;
        out[i] = (hi << 4) | lo;
    }
    Some(out)
}

fn hex_nibble(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

/// Convert any type that implements `noodles_sam::alignment::Record` into
/// our [`AlignmentRecord`]. Filters out unmapped / secondary / supplementary
/// records per the trait contract.
pub(crate) fn convert_trait_record<R: noodles_sam::alignment::Record>(
    record: &R,
) -> Option<AlignmentRecord> {
    let flags = record.flags().ok()?;
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let start = record.alignment_start()?.ok()?.get() as u32;
    let end = record.alignment_end()?.ok()?.get() as u32;
    let mapq = record
        .mapping_quality()
        .and_then(|r| r.ok())
        .map(u8::from)
        .unwrap_or(0);
    Some(AlignmentRecord {
        start,
        end,
        mapq,
        is_primary: true,
    })
}

/// Convert our [`Region`] to a `noodles_core::Region`. Shared by
/// `FileAlignmentSource`, `CursorBamSource`, and any future source that
/// queries noodles readers by region.
pub(crate) fn region_to_noodles(
    region: &crate::reader::Region,
) -> crate::Result<noodles_core::Region> {
    let start = noodles_core::Position::try_from(region.start as usize)
        .map_err(|_| crate::Error::Internal(format!("invalid region start: {}", region.start)))?;
    let end = noodles_core::Position::try_from(region.end as usize)
        .map_err(|_| crate::Error::Internal(format!("invalid region end: {}", region.end)))?;
    Ok(noodles_core::Region::new(region.chrom.clone(), start..=end))
}
