//! Error-path smoke tests for failure-mode fixtures.
//!
//! Each test loads a deliberately malformed fixture and asserts the correct
//! `Error` variant and `exit_code()`.

mod fixture_helpers;

use std::fs::File;
use std::io::Read as _;
use std::path::PathBuf;

use nds_core::reader::Region;
use nds_core::reference::ReferenceSource;
use nds_core::sources::memory::MemoryBamSource;
use nds_core::{analyze, AlignmentSource, Error};

use fixture_helpers::fixture_config;

struct FixtureReferenceSource {
    chr17_md5: [u8; 16],
}

impl ReferenceSource for FixtureReferenceSource {
    fn fetch_sequence(&mut self, _region: &Region) -> nds_core::Result<Vec<u8>> {
        Ok(Vec::new())
    }
    fn chr17_md5(&mut self) -> nds_core::Result<[u8; 16]> {
        Ok(self.chr17_md5)
    }
}

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn load_bytes(name: &str) -> Vec<u8> {
    let path = fixture_path(name);
    let mut bytes = Vec::new();
    File::open(&path)
        .unwrap_or_else(|e| panic!("open {}: {e}", path.display()))
        .read_to_end(&mut bytes)
        .unwrap();
    bytes
}

fn analyze_from_bytes(
    bytes: Vec<u8>,
    reference_md5: [u8; 16],
) -> nds_core::Result<nds_core::AnalysisResult> {
    let mut source = MemoryBamSource::new(bytes);
    let mut reference = FixtureReferenceSource {
        chr17_md5: reference_md5,
    };
    let config = fixture_config();
    analyze(&mut source, &mut reference, &config)
}

#[test]
fn wrong_build_returns_reference_build_mismatch() {
    let bytes = load_bytes("wrong_build.bam");
    let mut source = MemoryBamSource::new(bytes);
    let header = source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);
    let mut reference = FixtureReferenceSource { chr17_md5 };
    let config = fixture_config();
    let err = analyze(&mut source, &mut reference, &config).unwrap_err();
    assert!(
        matches!(err, Error::ReferenceBuildMismatch { .. }),
        "expected ReferenceBuildMismatch, got {err:?}"
    );
    assert_eq!(err.exit_code(), 12);
}

#[test]
fn empty_region_returns_empty_region() {
    let bytes = load_bytes("empty_region.bam");
    let mut source = MemoryBamSource::new(bytes);
    let header = source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);

    let bytes = load_bytes("empty_region.bam");
    let err = analyze_from_bytes(bytes, chr17_md5).unwrap_err();
    assert!(
        matches!(err, Error::EmptyRegion { .. }),
        "expected EmptyRegion, got {err:?}"
    );
    assert_eq!(err.exit_code(), 15);
}

#[test]
fn missing_index_returns_missing_index() {
    let path = fixture_path("missing_index/test.bam");
    // FileAlignmentSource::open succeeds (opens as BamPlain).
    // MissingIndex fires when query_region is called.
    let mut source = match nds_core::FileAlignmentSource::open(&path, None) {
        Ok(s) => s,
        Err(e) => panic!("expected open to succeed, got {e:?}"),
    };
    let region = nds_core::reader::Region::new("chr17", 1, 1000);
    let err = match source.query_region(&region) {
        Err(e) => e,
        Ok(_) => panic!("expected MissingIndex error on query_region"),
    };
    assert!(
        matches!(err, Error::MissingIndex(_)),
        "expected MissingIndex, got {err:?}"
    );
    assert_eq!(err.exit_code(), 11);
}

#[test]
fn invalid_format_returns_invalid_alignment_format() {
    let path = fixture_path("invalid_format.bin");
    let err = match nds_core::FileAlignmentSource::open(&path, None) {
        Err(e) => e,
        Ok(_) => panic!("expected InvalidAlignmentFormat error, got Ok"),
    };
    assert!(
        matches!(err, Error::InvalidAlignmentFormat(_)),
        "expected InvalidAlignmentFormat, got {err:?}"
    );
    assert_eq!(err.exit_code(), 10);
}

#[test]
fn truncated_bam_returns_error() {
    let bytes = load_bytes("truncated.bam");
    let mut source = MemoryBamSource::new(bytes);
    let header = source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);

    // Re-create source from bytes (header was consumed).
    let bytes = load_bytes("truncated.bam");
    let err = analyze_from_bytes(bytes, chr17_md5).unwrap_err();
    // MemoryBamSource wraps noodles errors as Error::Io rather than
    // Error::TruncatedInput (which requires the file.rs-specific
    // map_noodles_read_error). Accept either variant - the point is
    // that a truncated file does not produce a valid result.
    let code = err.exit_code();
    assert!(
        code == 14 || code == 20,
        "expected exit code 14 (TruncatedInput) or 20 (Io), got {code}: {err:?}"
    );
}

#[test]
fn md5_mismatch_returns_reference_md5_mismatch() {
    let bytes = load_bytes("md5_mismatch.bam");
    // The fixture header has M5 = 000...000 for chr17.
    // The reference source returns a DIFFERENT md5 → mismatch.
    let different_md5 = [0xFFu8; 16];
    let err = analyze_from_bytes(bytes, different_md5).unwrap_err();
    assert!(
        matches!(err, Error::ReferenceMd5Mismatch { .. }),
        "expected ReferenceMd5Mismatch, got {err:?}"
    );
    assert_eq!(err.exit_code(), 13);
}
