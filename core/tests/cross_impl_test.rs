//! Cross-implementation test.
//!
//! Compares FileAlignmentSource vs CursorBamSource on the same fixture
//! to prove the OffsetCursor path produces byte-identical results.

mod fixture_helpers;

use std::fs::File;
use std::io::Read as _;
use std::path::PathBuf;

use nds_core::reader::Region;
use nds_core::reference::ReferenceSource;
use nds_core::{analyze, AlignmentSource, CursorBamSource, FileAlignmentSource};

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
    File::open(&path).unwrap().read_to_end(&mut bytes).unwrap();
    bytes
}

#[test]
fn file_source_and_cursor_source_produce_identical_results() {
    let bam_name = "cn3_30x.bam";
    let bai_name = "cn3_30x.bam.bai";
    let config = fixture_config();

    // Path A: FileAlignmentSource (file-backed indexed reader).
    let path = fixture_path(bam_name);
    let mut file_source = FileAlignmentSource::open(&path, None).unwrap();
    let header = file_source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);
    let mut ref_a = FixtureReferenceSource { chr17_md5 };
    let result_a = analyze(&mut file_source, &mut ref_a, &config).unwrap();

    // Path B: CursorBamSource (OffsetCursor-backed indexed reader).
    let bam_bytes = load_bytes(bam_name);
    let bai_bytes = load_bytes(bai_name);
    let mut cursor_source = CursorBamSource::new(bam_bytes, bai_bytes).unwrap();
    let mut ref_b = FixtureReferenceSource { chr17_md5 };
    let result_b = analyze(&mut cursor_source, &mut ref_b, &config).unwrap();

    // AnalysisResult derives PartialEq - this is a field-level comparison
    // of depth ratio, boundary windows, and interpretation.
    assert_eq!(
        result_a, result_b,
        "FileAlignmentSource and CursorBamSource must produce identical AnalysisResults"
    );
}
