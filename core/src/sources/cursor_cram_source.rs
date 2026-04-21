//! `CursorCramSource` - an `AlignmentSource` wrapping an `OffsetCursor`
//! for indexed CRAM queries over prefetched byte ranges.
//!
//! Supports optional reference sequence bytes for full CRAM decoding.
//! When reference bytes are provided (browser path with bundled
//! sub-sequences), noodles-cram decodes records fully. When not provided,
//! sequence decoding may fail but alignment positions are still available
//! (depends on the CRAM's encoding - some files work without reference,
//! others don't).

use std::io::{BufReader, Cursor};

use crate::reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region};
use crate::sources::offset_cursor::OffsetCursor;
use crate::sources::{convert_trait_record, region_to_noodles, sam_header_to_alignment_header};
use crate::{Error, Result};

/// A pair of (chromosome_name, raw_sequence_bytes) for building the
/// noodles reference repository.
pub type ReferenceSequenceEntry = (String, Vec<u8>);

/// In-memory indexed CRAM source backed by an `OffsetCursor`.
pub struct CursorCramSource {
    reader: noodles_cram::io::Reader<BufReader<OffsetCursor>>,
    crai_index: noodles_cram::crai::Index,
    header: Option<noodles_sam::Header>,
    cached_header: Option<AlignmentHeader>,
}

impl CursorCramSource {
    /// Construct from a pre-built OffsetCursor + CRAI bytes + optional reference.
    pub fn from_cursor_with_ref(
        cursor: OffsetCursor,
        crai_bytes: Vec<u8>,
        reference_seqs: Vec<ReferenceSequenceEntry>,
    ) -> Result<Self> {
        let crai_index: noodles_cram::crai::Index =
            noodles_cram::crai::io::Reader::new(Cursor::new(crai_bytes))
                .read_index()
                .map_err(|e| {
                    Error::Internal(format!("failed to parse CRAI index from bytes: {e}"))
                })?;

        let mut builder = noodles_cram::io::reader::Builder::default();

        if !reference_seqs.is_empty() {
            let records: Vec<noodles_fasta::Record> = reference_seqs
                .into_iter()
                .map(|(name, seq_bytes)| {
                    let def = noodles_fasta::record::Definition::new(name, None);
                    let seq = noodles_fasta::record::Sequence::from(seq_bytes);
                    noodles_fasta::Record::new(def, seq)
                })
                .collect();
            let repo = noodles_fasta::repository::Repository::new(records);
            builder = builder.set_reference_sequence_repository(repo);
        }

        let reader = builder.build_from_reader(BufReader::new(cursor));
        Ok(Self {
            reader,
            crai_index,
            header: None,
            cached_header: None,
        })
    }

    /// Construct from a pre-built OffsetCursor + CRAI bytes (no reference).
    pub fn from_cursor(cursor: OffsetCursor, crai_bytes: Vec<u8>) -> Result<Self> {
        Self::from_cursor_with_ref(cursor, crai_bytes, Vec::new())
    }

    /// Construct from sparse byte ranges + CRAI bytes (no reference).
    pub fn from_ranges(ranges: Vec<(u64, Vec<u8>)>, crai_bytes: Vec<u8>) -> Result<Self> {
        Self::from_cursor(OffsetCursor::new(ranges), crai_bytes)
    }

    /// Construct from full contiguous CRAM + CRAI byte buffers (no reference).
    pub fn new(cram_bytes: Vec<u8>, crai_bytes: Vec<u8>) -> Result<Self> {
        Self::from_ranges(vec![(0, cram_bytes)], crai_bytes)
    }
}

impl AlignmentSource for CursorCramSource {
    fn read_header(&mut self) -> Result<AlignmentHeader> {
        if let Some(ref cached) = self.cached_header {
            return Ok(cached.clone());
        }
        let h = self
            .reader
            .read_header()
            .map_err(|e| Error::Io(std::io::Error::other(format!("read CRAM header: {e}"))))?;
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
            .query(header_ref, &self.crai_index, &noodles_region)
            .map_err(Error::Io)?;
        let iter = query.filter_map(|result| match result {
            Ok(record) => convert_trait_record(&record).map(Ok),
            Err(e) => Some(Err(Error::Io(std::io::Error::other(format!(
                "CRAM record decode: {e}"
            ))))),
        });
        Ok(Box::new(iter))
    }
}
