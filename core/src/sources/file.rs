//! File-backed [`AlignmentSource`] via noodles.
//!
//! Supports BAM (self-contained) and CRAM (requires a reference fasta).
//! Format auto-detected by sniffing the first 4 bytes.
//!
//! **Important:** gated `#[cfg(not(target_arch = "wasm32"))]` at the lib.rs
//! declaration; does not compile for the WASM build.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

use crate::reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region};
use crate::{Error, Result};

const CRAM_MAGIC: &[u8; 4] = b"CRAM";
const BGZF_MAGIC: &[u8; 2] = &[0x1f, 0x8b];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FileFormat {
    Bam,
    Cram,
}

fn detect_format(path: &Path) -> Result<FileFormat> {
    let mut file = File::open(path).map_err(Error::Io)?;
    let mut head = [0u8; 4];
    let n = file.read(&mut head).map_err(Error::Io)?;
    if n < 2 {
        return Err(Error::InvalidAlignmentFormat(path.display().to_string()));
    }
    if n >= 4 && &head == CRAM_MAGIC {
        return Ok(FileFormat::Cram);
    }
    if &head[..2] == BGZF_MAGIC {
        return Ok(FileFormat::Bam);
    }
    Err(Error::InvalidAlignmentFormat(path.display().to_string()))
}

/// Look for the companion index. Checks `<input>.<ext>` first
/// (e.g. `foo.bam.bai`), then `<input>` with the extension swapped
/// (e.g. `foo.bai`). Returns the first existing path.
fn detect_index(input: &Path, ext: &str) -> Option<PathBuf> {
    let appended = {
        let mut s = input.as_os_str().to_os_string();
        s.push(".");
        s.push(ext);
        PathBuf::from(s)
    };
    if appended.exists() {
        return Some(appended);
    }
    let swapped = input.with_extension(ext);
    if swapped.exists() {
        return Some(swapped);
    }
    None
}

/// Internal reader variants. Each variant caches its parsed header so
/// `read_header` can be called repeatedly without re-reading the file.
enum Inner {
    /// BAM with a `.bai` index. Supports indexed `query_region`.
    BamIndexed {
        reader: noodles_bam::io::IndexedReader<noodles_bgzf::io::Reader<BufReader<File>>>,
        header: Option<noodles_sam::Header>,
        cached_header: Option<AlignmentHeader>,
    },
    /// BAM without a `.bai` index. `read_header` works; `query_region`
    /// returns `Error::MissingIndex`. The inner reader is bgzf-wrapped
    /// because `noodles_bam::io::Reader::new` auto-wraps its argument in a
    /// BGZF reader.
    BamPlain {
        reader: noodles_bam::io::Reader<noodles_bgzf::io::Reader<BufReader<File>>>,
        header: Option<noodles_sam::Header>,
        cached_header: Option<AlignmentHeader>,
    },
    /// CRAM with or without an index. `query_region` loads the index lazily
    /// on first call and caches it in `cached_index`; if `index_path` is
    /// `None` it returns `Error::MissingIndex` immediately.
    Cram {
        reader: noodles_cram::io::Reader<BufReader<File>>,
        #[allow(dead_code)] // held so the repository outlives the reader
        repository: noodles_fasta::repository::Repository,
        index_path: Option<PathBuf>,
        cached_index: Option<noodles_cram::crai::Index>,
        header: Option<noodles_sam::Header>,
        cached_header: Option<AlignmentHeader>,
    },
}

/// File-backed alignment source.
pub struct FileAlignmentSource {
    inner: Inner,
}

impl FileAlignmentSource {
    /// Open an alignment file. `reference` is required for CRAM and ignored
    /// for BAM.
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self> {
        let format = detect_format(path)?;
        match format {
            FileFormat::Bam => Self::open_bam(path),
            FileFormat::Cram => Self::open_cram(path, reference),
        }
    }

    fn open_bam(path: &Path) -> Result<Self> {
        if let Some(index_path) = detect_index(path, "bai") {
            let index = noodles_bam::bai::fs::read(&index_path).map_err(|e| match e.kind() {
                std::io::ErrorKind::NotFound => {
                    Error::MissingIndex(index_path.display().to_string())
                }
                _ => Error::Io(e),
            })?;
            let file = File::open(path).map_err(Error::Io)?;
            let reader = noodles_bam::io::indexed_reader::Builder::default()
                .set_index(index)
                .build_from_reader(BufReader::new(file))
                .map_err(Error::Io)?;
            Ok(Self {
                inner: Inner::BamIndexed {
                    reader,
                    header: None,
                    cached_header: None,
                },
            })
        } else {
            let file = File::open(path).map_err(Error::Io)?;
            let reader = noodles_bam::io::Reader::new(BufReader::new(file));
            Ok(Self {
                inner: Inner::BamPlain {
                    reader,
                    header: None,
                    cached_header: None,
                },
            })
        }
    }

    fn open_cram(path: &Path, reference: Option<&Path>) -> Result<Self> {
        let reference = reference.ok_or_else(|| {
            Error::InvalidAlignmentFormat(format!(
                "{}: CRAM input requires a reference fasta (pass --reference or set REF_PATH)",
                path.display()
            ))
        })?;

        let fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(reference)
            .map_err(Error::Io)?;
        let adapter = noodles_fasta::repository::adapters::IndexedReader::new(fasta_reader);
        let repository = noodles_fasta::repository::Repository::new(adapter);

        let file = File::open(path).map_err(Error::Io)?;
        let reader = noodles_cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_reader(BufReader::new(file));
        let index_path = detect_index(path, "crai");

        Ok(Self {
            inner: Inner::Cram {
                reader,
                repository,
                index_path,
                cached_index: None,
                header: None,
                cached_header: None,
            },
        })
    }
}

// Shared SAM-header and record conversion helpers live in `sources/mod.rs`
// so both `FileAlignmentSource` and `MemoryBamSource` can use them without
// `memory.rs` having to pull in the file module (which is cfg-gated off on
// wasm32). Re-exported through `use` so the tests below can call them via
// `super::hex_to_md5_bytes` exactly as they did before the refactor.
#[cfg(test)]
pub(super) use crate::sources::hex_to_md5_bytes;
pub(super) use crate::sources::{
    convert_trait_record, region_to_noodles, sam_header_to_alignment_header,
};

fn map_noodles_read_error(e: std::io::Error) -> Error {
    match e.kind() {
        std::io::ErrorKind::UnexpectedEof => Error::TruncatedInput,
        _ => Error::Io(e),
    }
}

impl AlignmentSource for FileAlignmentSource {
    fn read_header(&mut self) -> Result<AlignmentHeader> {
        match &mut self.inner {
            Inner::BamIndexed {
                reader,
                header,
                cached_header,
            } => {
                if let Some(ref cached) = cached_header {
                    return Ok(cached.clone());
                }
                let h = reader.read_header().map_err(map_noodles_read_error)?;
                let converted = sam_header_to_alignment_header(&h);
                *header = Some(h);
                *cached_header = Some(converted.clone());
                Ok(converted)
            }
            Inner::BamPlain {
                reader,
                header,
                cached_header,
            } => {
                if let Some(ref cached) = cached_header {
                    return Ok(cached.clone());
                }
                let h = reader.read_header().map_err(map_noodles_read_error)?;
                let converted = sam_header_to_alignment_header(&h);
                *header = Some(h);
                *cached_header = Some(converted.clone());
                Ok(converted)
            }
            Inner::Cram {
                reader,
                header,
                cached_header,
                ..
            } => {
                if let Some(ref cached) = cached_header {
                    return Ok(cached.clone());
                }
                let h = reader.read_header().map_err(map_noodles_read_error)?;
                let converted = sam_header_to_alignment_header(&h);
                *header = Some(h);
                *cached_header = Some(converted.clone());
                Ok(converted)
            }
        }
    }

    #[allow(clippy::type_complexity)]
    fn query_region<'a>(
        &'a mut self,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<AlignmentRecord>> + 'a>> {
        self.read_header()?;
        let noodles_region = region_to_noodles(region)?;

        match &mut self.inner {
            Inner::BamIndexed { reader, header, .. } => {
                let header_ref = header.as_ref().expect("header cached above");
                let query = reader
                    .query(header_ref, &noodles_region)
                    .map_err(Error::Io)?;
                // BAM Query does NOT implement Iterator directly; `.records()`
                // converts it into one. (CRAM Query does implement Iterator
                // and doesn't need this call.)
                let iter = query.records().filter_map(|result| match result {
                    Ok(record) => convert_trait_record(&record).map(Ok),
                    Err(e) => Some(Err(map_noodles_read_error(e))),
                });
                Ok(Box::new(iter))
            }
            Inner::BamPlain { .. } => Err(Error::MissingIndex(
                "no .bai next to the BAM input (expected <input>.bai or <input.bam>.bai)".into(),
            )),
            Inner::Cram {
                reader,
                index_path,
                cached_index,
                header,
                ..
            } => {
                // Lazily load and cache the CRAI on first call. The cached
                // index lives with the reader so the query iterator's borrow
                // against it doesn't escape the function.
                if cached_index.is_none() {
                    let path = index_path.clone().ok_or_else(|| {
                        Error::MissingIndex(
                            "no .crai next to the CRAM input (expected <input>.crai or <input.cram>.crai)".into(),
                        )
                    })?;
                    let index =
                        noodles_cram::crai::fs::read(&path).map_err(|e| match e.kind() {
                            std::io::ErrorKind::NotFound => {
                                Error::MissingIndex(path.display().to_string())
                            }
                            _ => Error::Io(e),
                        })?;
                    *cached_index = Some(index);
                }
                let index = cached_index.as_ref().expect("just cached");
                let header_ref = header.as_ref().expect("header cached above");
                let query = reader
                    .query(header_ref, index, &noodles_region)
                    .map_err(Error::Io)?;
                let iter = query.filter_map(|result| match result {
                    Ok(record) => convert_trait_record(&record).map(Ok),
                    Err(e) => Some(Err(map_noodles_read_error(e))),
                });
                Ok(Box::new(iter))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_minimal_bam(path: &Path) -> std::io::Result<[u8; 16]> {
        use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
        use std::num::NonZeroUsize;

        const MD5_HEX: &str = "abcdef0123456789abcdef0123456789";
        let md5_bytes = super::hex_to_md5_bytes(MD5_HEX.as_bytes()).unwrap();

        let mut ref_seq: Map<ReferenceSequence> =
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(83_257_441usize).unwrap());
        let tag = noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
        ref_seq.other_fields_mut().insert(tag, MD5_HEX.into());

        let header = noodles_sam::Header::builder()
            .add_reference_sequence("chr17", ref_seq)
            .build();

        let file = File::create(path)?;
        let mut writer = noodles_bam::io::Writer::new(file);
        writer.write_header(&header)?;
        Ok(md5_bytes)
    }

    #[test]
    fn file_alignment_source_round_trips_a_minimal_bam() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("spike.bam");
        let expected_md5 = write_minimal_bam(&bam_path).unwrap();

        let mut source = FileAlignmentSource::open(&bam_path, None).unwrap();
        let header = source.read_header().unwrap();
        assert_eq!(header.sequences.len(), 1);
        let seq = &header.sequences[0];
        assert_eq!(seq.name, "chr17");
        assert_eq!(seq.length, 83_257_441);
        assert_eq!(seq.md5, Some(expected_md5));
    }

    fn open_expect_err(path: &Path, reference: Option<&Path>) -> Error {
        match FileAlignmentSource::open(path, reference) {
            Ok(_) => panic!("expected FileAlignmentSource::open to fail"),
            Err(e) => e,
        }
    }

    #[test]
    fn file_alignment_source_rejects_invalid_format() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("garbage.bam");
        {
            let mut f = File::create(&path).unwrap();
            f.write_all(&[0u8; 16]).unwrap();
        }
        let err = open_expect_err(&path, None);
        assert!(
            matches!(err, Error::InvalidAlignmentFormat(_)),
            "expected InvalidAlignmentFormat, got {err:?}"
        );
    }

    #[test]
    fn file_alignment_source_rejects_cram_without_reference() {
        let tmp = tempfile::tempdir().unwrap();
        let path = tmp.path().join("fake.cram");
        {
            let mut f = File::create(&path).unwrap();
            f.write_all(b"CRAM\x03").unwrap();
        }
        let err = open_expect_err(&path, None);
        let msg = err.to_string();
        assert!(
            msg.contains("reference"),
            "error message should mention reference: {msg}"
        );
    }

    #[test]
    fn file_alignment_source_unindexed_bam_query_returns_missing_index() {
        let tmp = tempfile::tempdir().unwrap();
        let bam_path = tmp.path().join("unindexed.bam");
        write_minimal_bam(&bam_path).unwrap();
        // No .bai index written next to it.
        let mut source = FileAlignmentSource::open(&bam_path, None).unwrap();
        let region = Region::new("chr17", 1000, 2000);
        let err = match source.query_region(&region) {
            Ok(_) => panic!("expected MissingIndex error"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::MissingIndex(_)),
            "expected MissingIndex, got {err:?}"
        );
        assert_eq!(err.exit_code(), 11);
    }

    #[test]
    fn hex_to_md5_bytes_rejects_bad_input() {
        assert_eq!(hex_to_md5_bytes(b""), None);
        assert_eq!(hex_to_md5_bytes(b"abc"), None);
        assert_eq!(hex_to_md5_bytes(b"abcdef0123456789abcdef012345678Z"), None);
        assert_eq!(
            hex_to_md5_bytes(b"abcdef0123456789abcdef0123456789"),
            Some([
                0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45,
                0x67, 0x89,
            ])
        );
    }

    #[test]
    fn hex_to_md5_bytes_is_case_insensitive() {
        assert_eq!(
            hex_to_md5_bytes(b"ABCDEF0123456789ABCDEF0123456789"),
            Some([
                0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45,
                0x67, 0x89,
            ])
        );
    }
}
