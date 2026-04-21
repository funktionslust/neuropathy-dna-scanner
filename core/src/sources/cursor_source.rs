//! `CursorBamSource` - an `AlignmentSource` wrapping an `OffsetCursor`
//! in a noodles `IndexedReader` for indexed BAM queries over prefetched
//! byte ranges.
//!
//! This is the WASM-path equivalent of `FileAlignmentSource`'s
//! `Inner::BamIndexed`. The cursor holds the BAM bytes in memory; the
//! BAI index is parsed from in-memory bytes. noodles' `IndexedReader`
//! issues seeks and reads against the cursor exactly as it would against
//! a file.

use std::io::Cursor;

use crate::reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region};
use crate::sources::offset_cursor::OffsetCursor;
use crate::sources::{convert_trait_record, region_to_noodles, sam_header_to_alignment_header};
use crate::{Error, Result};

/// In-memory indexed BAM source backed by an `OffsetCursor`. Use
/// [`Self::new`] for a single contiguous buffer, [`Self::from_ranges`] for
/// sparse `(file_offset, bytes)` inputs, or [`Self::from_cursor`] to share
/// an existing [`OffsetCursor`].
pub struct CursorBamSource {
    reader:
        noodles_bam::io::IndexedReader<noodles_bgzf::io::Reader<std::io::BufReader<OffsetCursor>>>,
    header: Option<noodles_sam::Header>,
    cached_header: Option<AlignmentHeader>,
}

impl CursorBamSource {
    /// Construct from a pre-built OffsetCursor + BAI bytes.
    pub fn from_cursor(cursor: OffsetCursor, bai_bytes: Vec<u8>) -> Result<Self> {
        let index: noodles_bam::bai::Index =
            noodles_bam::bai::io::Reader::new(Cursor::new(bai_bytes))
                .read_index()
                .map_err(|e| {
                    Error::Internal(format!("failed to parse BAI index from bytes: {e}"))
                })?;
        let reader = noodles_bam::io::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_reader(std::io::BufReader::new(cursor))
            .map_err(Error::Io)?;
        Ok(Self {
            reader,
            header: None,
            cached_header: None,
        })
    }

    /// Construct from sparse byte ranges fetched via `File.slice()`.
    /// Each range is `(file_offset, bytes)`.
    pub fn from_ranges(ranges: Vec<(u64, Vec<u8>)>, bai_bytes: Vec<u8>) -> Result<Self> {
        Self::from_cursor(OffsetCursor::new(ranges), bai_bytes)
    }

    /// Construct from full contiguous BAM + BAI byte buffers.
    pub fn new(bam_bytes: Vec<u8>, bai_bytes: Vec<u8>) -> Result<Self> {
        Self::from_ranges(vec![(0, bam_bytes)], bai_bytes)
    }
}

impl AlignmentSource for CursorBamSource {
    fn read_header(&mut self) -> Result<AlignmentHeader> {
        if let Some(ref cached) = self.cached_header {
            return Ok(cached.clone());
        }
        let h = self
            .reader
            .read_header()
            .map_err(|e| Error::Io(std::io::Error::other(format!("read BAM header: {e}"))))?;
        let converted = sam_header_to_alignment_header(&h);
        self.header = Some(h);
        self.cached_header = Some(converted.clone());
        Ok(converted)
    }

    fn query_region<'a>(
        &'a mut self,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<AlignmentRecord>> + 'a>> {
        self.read_header()?;
        let noodles_region = region_to_noodles(region)?;
        let header_ref = self.header.as_ref().expect("header cached above");
        let query = self
            .reader
            .query(header_ref, &noodles_region)
            .map_err(Error::Io)?;
        let iter = query.records().filter_map(|result| match result {
            Ok(record) => convert_trait_record(&record).map(Ok),
            Err(e) => Some(Err(Error::Io(std::io::Error::other(format!(
                "BAM record decode: {e}"
            ))))),
        });
        Ok(Box::new(iter))
    }
}
