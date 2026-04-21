//! Smoke tests for the baseline BAM fixtures.
//!
//! Each test loads a committed fixture via `MemoryBamSource`, runs
//! `analyze()` with `fixture_config()`, and asserts the expected
//! CN + confidence + subtype.

mod fixture_helpers;

use std::fs::File;
use std::io::Read as _;
use std::path::PathBuf;

use nds_core::interpret::{Confidence, Subtype};
use nds_core::reader::Region;
use nds_core::reference::ReferenceSource;
use nds_core::sources::memory::MemoryBamSource;
use nds_core::{analyze, AlignmentSource};

use fixture_helpers::fixture_config;

/// Inline reference source for integration tests - returns the BAM
/// header's own chr17 M5 tag, same pattern as `real_data_probe`.
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

fn load_and_analyze(bam_name: &str) -> nds_core::Result<nds_core::AnalysisResult> {
    let path = fixture_path(bam_name);
    let mut bytes = Vec::new();
    File::open(&path)
        .unwrap_or_else(|e| {
            panic!(
                "open {}: {e} - run `cargo xtask synth-fixtures` first",
                path.display()
            )
        })
        .read_to_end(&mut bytes)
        .unwrap();

    let mut source = MemoryBamSource::new(bytes);
    let header = source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);
    let mut reference = FixtureReferenceSource { chr17_md5 };
    let config = fixture_config();

    analyze(&mut source, &mut reference, &config)
}

#[test]
fn cn2_fixture_produces_normal_result() {
    let result = load_and_analyze("cn2_30x.bam").unwrap();
    assert_eq!(result.depth.estimated_cn, 2, "expected CN=2");
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(result.interpretation.subtype_call, None);
    assert!(result.interpretation.atypical_flags.is_empty());
}

#[test]
fn cn3_fixture_produces_cmt1a_result() {
    let result = load_and_analyze("cn3_30x.bam").unwrap();
    assert_eq!(result.depth.estimated_cn, 3, "expected CN=3");
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(
        result.interpretation.subtype_call,
        Some(Subtype::CMT1A),
        "expected Subtype::CMT1A"
    );
}

#[test]
fn cn1_fixture_produces_hnpp_result() {
    let result = load_and_analyze("cn1_30x.bam").unwrap();
    assert_eq!(result.depth.estimated_cn, 1, "expected CN=1");
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(
        result.interpretation.subtype_call,
        Some(Subtype::HNPP),
        "expected Subtype::HNPP"
    );
}

// --- Coverage-floor fixtures ---

#[test]
fn cn3_low_20x_fixture_produces_low_confidence() {
    let result = load_and_analyze("cn3_low_20x.bam").unwrap();
    assert_eq!(result.depth.estimated_cn, 3, "expected CN=3");
    assert_eq!(
        result.interpretation.confidence,
        Confidence::Low,
        "expected Confidence::Low at 20x autosomal"
    );
    assert!(
        result.interpretation.copy_number.is_some(),
        "low-confidence should still return a copy number"
    );
}

#[test]
fn cn3_floor_10x_fixture_produces_refused() {
    let result = load_and_analyze("cn3_floor_10x.bam").unwrap();
    assert_eq!(
        result.interpretation.confidence,
        Confidence::Refused,
        "expected Confidence::Refused at 10x autosomal"
    );
    assert_eq!(
        result.interpretation.copy_number, None,
        "refused should not return a copy number"
    );
}

#[test]
fn cn2_26x_fixture_produces_full_confidence() {
    let result = load_and_analyze("cn2_26x.bam").unwrap();
    assert_eq!(result.depth.estimated_cn, 2, "expected CN=2");
    assert_eq!(
        result.interpretation.confidence,
        Confidence::Full,
        "expected Confidence::Full at 26x autosomal (above 25x threshold with margin for jitter)"
    );
}
