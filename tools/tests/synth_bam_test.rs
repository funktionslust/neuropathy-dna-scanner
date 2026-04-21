//! Round-trip integration test for synth_bam: generate a CN=2 BAM,
//! read it back via MemoryBamSource, run depth_ratio(), and assert the
//! ratio is approximately 1.0 (CN=2).

use std::fs::File;
use std::io::Read as _;

use nds_core::depth::depth_ratio;
use nds_core::reader::Region;
use nds_core::sources::memory::MemoryBamSource;
use nds_core::AnalysisConfig;

#[test]
fn synth_bam_cn2_round_trip_via_memory_source() {
    // Use a small chr2 control region (500 kb instead of the default 20 MB)
    // so the test generates ~100k reads instead of ~4M. MemoryBamSource does
    // naive full scans per query, so read count dominates runtime.
    let config = AnalysisConfig {
        chr2_control_region: Region::new("chr2", 50_000_000, 50_500_000),
        ..AnalysisConfig::default()
    };

    let specs = vec![
        nds_tools::synth::CoverageSpec {
            chrom: config.pmp22_region.chrom.clone(),
            start: config.pmp22_region.start,
            end: config.pmp22_region.end,
            depth: 30,
        },
        nds_tools::synth::CoverageSpec {
            chrom: config.chr2_control_region.chrom.clone(),
            start: config.chr2_control_region.start,
            end: config.chr2_control_region.end,
            depth: 30,
        },
    ];

    let tmp = tempfile::tempdir().unwrap();
    let bam_path = tmp.path().join("cn2_roundtrip.bam");

    // Generate the BAM (we skip indexing since MemoryBamSource doesn't need it).
    let header = nds_tools::synth::build_grch38_header();
    let reads = nds_tools::synth::generate_reads(&header, &specs, 42);
    nds_tools::synth::write_bam(&bam_path, &header, &reads).unwrap();

    // Slurp bytes and analyze via MemoryBamSource.
    let mut bytes = Vec::new();
    File::open(&bam_path)
        .unwrap()
        .read_to_end(&mut bytes)
        .unwrap();

    let mut source = MemoryBamSource::new(bytes);
    let result = depth_ratio(&mut source, &config).unwrap();

    assert!(
        (result.ratio - 1.0).abs() < 0.05,
        "expected ratio ~1.0 (CN=2), got {:.3}",
        result.ratio
    );
    assert_eq!(
        result.estimated_cn, 2,
        "expected CN=2, got {}",
        result.estimated_cn
    );
}
