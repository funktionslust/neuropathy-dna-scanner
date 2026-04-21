//! Integration tests exercising the FileAlignmentSource → analyze() path
//! on committed BAM fixtures.
//!
//! These prove the indexed BAM reader, header parsing, and region query
//! pipeline works end-to-end. Algorithm assertions (CN, confidence,
//! subtype) are the same as the MemoryBamSource smoke tests  - the value
//! here is testing the I/O layer, not the algorithm.

mod fixture_helpers;

use nds_core::interpret::{Confidence, Subtype};
use nds_core::reader::Region;
use nds_core::reference::ReferenceSource;
use nds_core::{analyze, AlignmentSource, FileAlignmentSource};
use std::path::PathBuf;

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

fn open_and_analyze(bam_name: &str) -> nds_core::AnalysisResult {
    let path = fixture_path(bam_name);
    let mut source = FileAlignmentSource::open(&path, None)
        .unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let header = source.read_header().unwrap();
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);
    let mut reference = FixtureReferenceSource { chr17_md5 };
    let config = fixture_config();
    analyze(&mut source, &mut reference, &config).unwrap()
}

#[test]
fn file_source_cn2_produces_normal() {
    let result = open_and_analyze("cn2_30x.bam");
    assert_eq!(result.depth.estimated_cn, 2);
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(result.interpretation.subtype_call, None);
}

#[test]
fn file_source_cn3_produces_cmt1a() {
    let result = open_and_analyze("cn3_30x.bam");
    assert_eq!(result.depth.estimated_cn, 3);
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(result.interpretation.subtype_call, Some(Subtype::CMT1A));
}

#[test]
fn file_source_cn1_produces_hnpp() {
    let result = open_and_analyze("cn1_30x.bam");
    assert_eq!(result.depth.estimated_cn, 1);
    assert_eq!(result.interpretation.confidence, Confidence::Full);
    assert_eq!(result.interpretation.subtype_call, Some(Subtype::HNPP));
}
