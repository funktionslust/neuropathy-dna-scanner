//! neuropathy-dna-scanner core library.
//!
//! Depth-ratio and boundary-scan algorithms over CRAM/BAM alignments, plus the
//! supported-subtypes matrix that is the single source of truth.
//! Target-agnostic: the same code runs in the native CLI binary and inside the
//! `wasm32-unknown-unknown` browser build.

pub mod boundary;
pub mod catalog;
pub mod depth;
pub mod error;
pub mod interpret;
pub mod matrix;
pub mod pileup;
pub mod pileup_cram;
pub mod reader;
pub mod reference;
pub mod reference_data;
pub mod report;

pub mod sources;

#[cfg(test)]
pub(crate) mod mock;

pub use boundary::{
    boundary_scan, BoundaryDirection, BoundaryRun, BoundaryScanResult, WindowDepth,
};
pub use depth::{depth_ratio, DepthRatio};
pub use error::{Error, Result};
pub use interpret::{interpret, AtypicalFlag, Confidence, Interpretation, Subtype};
pub use matrix::{
    render_html, render_markdown, render_text, MatrixEntry, SubtypeStatus, DISCLAIMER, SUPPORTED,
};
pub use reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region, SequenceInfo};
pub use reference::ReferenceSource;
pub use report::{render_report, ReportMeta};

pub use sources::cursor_cram_source::CursorCramSource;
pub use sources::cursor_source::CursorBamSource;
#[cfg(not(target_arch = "wasm32"))]
pub use sources::file::FileAlignmentSource;
pub use sources::memory::MemoryBamSource;

/// Configuration for a single analysis run.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisConfig {
    /// PMP22 gene region for the depth-ratio numerator. Defaults to the
    /// MANE-select transcript coordinates on GRCh38.
    pub pmp22_region: Region,
    /// Chromosome 2 control region for the depth-ratio denominator.
    pub chr2_control_region: Region,
    /// Region scanned by `boundary_scan` in 10 kb windows. Extends to
    /// `chr17:21,000,000` so the SMS-REP proximal element is covered for
    /// PTLS/SMS/YUHAL detection.
    pub boundary_scan_region: Region,
    /// Autosomal coverage below which no CN call is returned.
    pub coverage_floor_refuse: f64,
    /// Autosomal coverage below which the call is marked low-confidence.
    pub coverage_floor_low: f64,
    /// Minimum boundary-run length (bases) for a run to be classified as the
    /// classical 1.4 Mb CMT1A duplication or HNPP deletion.
    pub min_classical_run_length: u32,
    /// Maximum boundary-run length (bases) for a run to be classified as
    /// classical.
    pub max_classical_run_length: u32,
    /// Maximum number of consecutive below-threshold windows allowed inside
    /// a boundary-scan run before the run is closed. Lets the run detector
    /// bridge CMT1A-REP segdup dips without merging distinct CNV events.
    /// Default: 5 (50 kb at the 10 kb window size), empirically validated
    /// on a CMT1A-positive CRAM.
    pub max_run_gap_windows: usize,
    /// RAI1 gene region for the PTLS/SMS depth-ratio check. Defaults to
    /// `chr17:17,681,458-17,811,453` (GRCh38, NCBI Gene 10743).
    pub rai1_region: Region,
    /// Normalised depth above which a window is classified as duplicated.
    /// Default: 1.3 (CN=3 has ratio ~1.5; 1.3 sits between CN=2 baseline
    /// of 1.0 and the expected CN=3 signal).
    pub duplication_threshold: f64,
    /// Normalised depth below which a window is classified as deleted.
    /// Default: 0.7 (CN=1 has ratio ~0.5; 0.7 sits between baseline 1.0
    /// and the expected CN=1 signal).
    pub deletion_threshold: f64,
    /// Tolerance in bases for matching a boundary run endpoint to a known
    /// REP element. Accounts for 10 kb window resolution and
    /// segmental-duplication-edge mapping noise.
    pub rep_match_tolerance: u32,
    /// Boundary run shorter than this is flagged as atypical rather than
    /// classified as a standard syndrome. Default: 500,000 (500 kb).
    /// Standard CMT1A is ~1.4 Mb; anything under 500 kb is too small for
    /// a recurrent NAHR event in this region.
    pub atypical_short_threshold: u32,
    /// Minimum fraction of reads with MAPQ >= 20 for a window to
    /// contribute to run detection. Windows below this threshold are
    /// paralog-rich zones where short-read depth is unreliable.
    /// Default: 0.5 (at least half the reads must map uniquely).
    pub min_unique_mapper_fraction: f64,
    /// Minimum run length in bases for a run to be reported. Runs
    /// shorter than this are discarded as depth noise or segdup
    /// artifacts. Default: 100,000 (100 kb). No clinically meaningful
    /// CMT-region CNV is shorter than this.
    pub min_reportable_run_length: u32,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            // NCBI Gene 5376, GRCh38.p14. Matches reference_data::grch38().genes[0].
            pmp22_region: Region::new("chr17", 15_229_779, 15_265_326),
            chr2_control_region: Region::new("chr2", 50_000_000, 70_000_000),
            boundary_scan_region: Region::new("chr17", 13_500_000, 21_000_000),
            coverage_floor_refuse: 15.0,
            coverage_floor_low: 25.0,
            min_classical_run_length: 1_000_000,
            max_classical_run_length: 2_000_000,
            max_run_gap_windows: 5,
            rai1_region: Region::new("chr17", 17_681_458, 17_811_453),
            duplication_threshold: 1.3,
            deletion_threshold: 0.7,
            rep_match_tolerance: 100_000,
            atypical_short_threshold: 500_000,
            min_unique_mapper_fraction: 0.5,
            min_reportable_run_length: 100_000,
        }
    }
}

impl AnalysisConfig {
    /// Validate config values for internal consistency. Returns an error
    /// describing the first invalid field found.
    pub fn validate(&self) -> Result<()> {
        if self.coverage_floor_refuse <= 0.0 {
            return Err(Error::InvalidConfig(format!(
                "coverage_floor_refuse must be positive, got {}",
                self.coverage_floor_refuse
            )));
        }
        if self.coverage_floor_low <= self.coverage_floor_refuse {
            return Err(Error::InvalidConfig(format!(
                "coverage_floor_low ({}) must exceed coverage_floor_refuse ({})",
                self.coverage_floor_low, self.coverage_floor_refuse
            )));
        }
        if self.duplication_threshold <= 1.0 {
            return Err(Error::InvalidConfig(format!(
                "duplication_threshold must be > 1.0, got {}",
                self.duplication_threshold
            )));
        }
        if self.deletion_threshold >= 1.0 || self.deletion_threshold <= 0.0 {
            return Err(Error::InvalidConfig(format!(
                "deletion_threshold must be between 0.0 and 1.0, got {}",
                self.deletion_threshold
            )));
        }
        if self.min_classical_run_length >= self.max_classical_run_length {
            return Err(Error::InvalidConfig(format!(
                "min_classical_run_length ({}) must be < max_classical_run_length ({})",
                self.min_classical_run_length, self.max_classical_run_length
            )));
        }
        if self.pmp22_region.start >= self.pmp22_region.end {
            return Err(Error::InvalidConfig(
                "pmp22_region is inverted (start >= end)".to_string(),
            ));
        }
        if self.min_unique_mapper_fraction < 0.0 || self.min_unique_mapper_fraction > 1.0 {
            return Err(Error::InvalidConfig(format!(
                "min_unique_mapper_fraction must be in [0.0, 1.0], got {}",
                self.min_unique_mapper_fraction
            )));
        }
        Ok(())
    }
}

/// Detected biological sex from chrX/chrY depth ratios.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

/// Detect sex from chrX and chrY depth-proxy ratios relative to autosomal baseline.
///
/// Uses two independent signals that must agree:
/// - **chrX/autosomal**: males ~0.5 (hemizygous), females ~1.0 (diploid)
/// - **chrY/chrX**: males ~0.3-0.7, females <0.10 (noise only)
///
/// The chrX/autosomal signal is the primary discriminator (2x gap, robust
/// against CRAM compression bias). chrY/chrX is confirmatory.
/// If the signals disagree, returns Unknown.
pub fn detect_sex(chrx_depth: f64, chry_depth: f64, autosomal_depth: f64) -> Sex {
    if autosomal_depth < 1e-10 || chrx_depth < 1e-10 {
        return Sex::Unknown;
    }
    let chrx_to_auto = chrx_depth / autosomal_depth;
    let chry_to_chrx = chry_depth / chrx_depth;

    let chrx_says_male = chrx_to_auto < 0.65;
    let chrx_says_female = chrx_to_auto > 0.80;
    let chry_says_male = chry_to_chrx > 0.15;
    let chry_says_female = chry_to_chrx < 0.15;

    if chrx_says_male && chry_says_male {
        Sex::Male
    } else if chrx_says_female && chry_says_female {
        Sex::Female
    } else {
        Sex::Unknown
    }
}

/// Full output of a single analysis run. Produced by [`analyze`]; consumed
/// by [`render_report`].
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisResult {
    pub depth: depth::DepthRatio,
    pub boundary: boundary::BoundaryScanResult,
    pub interpretation: interpret::Interpretation,
}

/// GRCh38 chr17 length in bases. Source: NCBI Assembly GRCh38.p14.
pub const GRCH38_CHR17_LENGTH: u64 = 83_257_441;
/// GRCh37 chr17 length in bases. Source: NCBI Assembly GRCh37.p13.
pub const GRCH37_CHR17_LENGTH: u64 = 81_195_210;

/// Detected reference genome build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceBuild {
    GRCh37,
    GRCh38,
}

/// Progress events emitted by [`analyze_with_progress`] at phase boundaries
/// and during long record-walking loops.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum Progress {
    /// Header has been read and the chr17 build/MD5 check passed.
    HeaderVerified,
    /// Starting a `query_region` walk for a specific labelled region.
    WalkStart { phase: &'static str, region: String },
    /// Per-N-records heartbeat from inside a walk. `records_seen` is the
    /// cumulative count in the current walk; `last_position` is the 1-based
    /// reference coordinate of the most recent record's start.
    /// Emitted every [`PROGRESS_TICK_INTERVAL`] records.
    WalkTick {
        phase: &'static str,
        records_seen: u64,
        last_position: u32,
    },
    /// A walk completed. `records_seen` is the final count.
    WalkDone {
        phase: &'static str,
        records_seen: u64,
    },
    /// `depth_ratio` finished and produced a numeric result.
    DepthRatioDone {
        pmp22_mean_depth: f64,
        autosomal_mean_depth: f64,
        ratio: f64,
        estimated_cn: u8,
    },
    /// `boundary_scan` finished.
    BoundaryScanDone {
        windows: usize,
        longest_run_length: Option<u32>,
    },
    /// `interpret` finished.
    InterpretDone {
        confidence: interpret::Confidence,
        subtype_call: Option<interpret::Subtype>,
    },
}

/// How often (in records) [`analyze_with_progress`] emits a [`Progress::WalkTick`].
pub const PROGRESS_TICK_INTERVAL: u64 = 50_000;

/// Run the full analysis pipeline: header verification, depth ratio,
/// boundary scan, interpretation. Returns an aggregated [`AnalysisResult`].
///
/// Coverage below the floor returns `Ok(AnalysisResult)` with
/// `confidence == Confidence::Refused` and `copy_number == None` — it is not
/// an error.
pub fn analyze<A, R>(
    source: &mut A,
    reference: &mut R,
    config: &AnalysisConfig,
) -> Result<AnalysisResult>
where
    A: reader::AlignmentSource,
    R: reference::ReferenceSource,
{
    analyze_with_progress(source, reference, config, |_| {})
}

/// Same as [`analyze`] but emits [`Progress`] events to a caller-supplied
/// callback at phase boundaries and during long record-walking loops.
pub fn analyze_with_progress<A, R, F>(
    source: &mut A,
    reference: &mut R,
    config: &AnalysisConfig,
    mut on_progress: F,
) -> Result<AnalysisResult>
where
    A: reader::AlignmentSource,
    R: reference::ReferenceSource,
    F: FnMut(Progress),
{
    config.validate()?;

    let header = source.read_header()?;
    let build = verify_build_and_md5(&header, reference)?;
    on_progress(Progress::HeaderVerified);

    let depth = depth::depth_ratio_with_progress(source, config, &mut on_progress)?;
    on_progress(Progress::DepthRatioDone {
        pmp22_mean_depth: depth.pmp22_mean_depth,
        autosomal_mean_depth: depth.autosomal_mean_depth,
        ratio: depth.ratio,
        estimated_cn: depth.estimated_cn,
    });

    let boundary = boundary::boundary_scan_with_progress(
        source,
        depth.autosomal_mean_depth,
        config,
        &mut on_progress,
    )?;
    on_progress(Progress::BoundaryScanDone {
        windows: boundary.windows.len(),
        longest_run_length: boundary.longest_run().map(|r| r.length),
    });

    let interpretation = interpret::interpret(&depth, &boundary, config, build);
    on_progress(Progress::InterpretDone {
        confidence: interpretation.confidence,
        subtype_call: interpretation.subtype_call,
    });

    Ok(AnalysisResult {
        depth,
        boundary,
        interpretation,
    })
}

/// Header reference-build and MD5 verification.
///
/// 1. chr17 must be present (case-insensitive `chr17` or `17`) with a
///    known assembly length; otherwise [`Error::ReferenceBuildMismatch`].
/// 2. If the header's chr17 has an `M5` tag, it must equal
///    `reference.chr17_md5()`; otherwise [`Error::ReferenceMd5Mismatch`].
///    Missing `M5` is not an error — real-world files frequently omit it.
fn verify_build_and_md5<R>(
    header: &reader::AlignmentHeader,
    reference: &mut R,
) -> Result<ReferenceBuild>
where
    R: reference::ReferenceSource,
{
    if header.sequences.is_empty() {
        return Err(Error::ReferenceBuildMismatch {
            got: "alignment header has no reference sequences (empty header)".to_string(),
        });
    }

    let chr17_seq = header
        .sequences
        .iter()
        .find(|s| s.name.eq_ignore_ascii_case("chr17") || s.name.eq_ignore_ascii_case("17"))
        .ok_or_else(|| Error::ReferenceBuildMismatch {
            got: "no chr17 or 17 sequence in the alignment header".to_string(),
        })?;

    let build = if chr17_seq.length == GRCH38_CHR17_LENGTH {
        ReferenceBuild::GRCh38
    } else if chr17_seq.length == GRCH37_CHR17_LENGTH {
        ReferenceBuild::GRCh37
    } else {
        return Err(Error::ReferenceBuildMismatch {
            got: format!(
                "chr17 length {} (expected GRCh38={} or GRCh37={})",
                chr17_seq.length, GRCH38_CHR17_LENGTH, GRCH37_CHR17_LENGTH
            ),
        });
    };

    if let Some(header_md5) = chr17_seq.md5 {
        let ref_md5 = reference.chr17_md5()?;
        if header_md5 != ref_md5 {
            return Err(Error::ReferenceMd5Mismatch {
                expected: md5_hex(header_md5),
                got: md5_hex(ref_md5),
            });
        }
    }

    Ok(build)
}

/// Convert an MD5 byte array to its 32-character lowercase hex string.
fn md5_hex(md5: [u8; 16]) -> String {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(32);
    for byte in md5 {
        write!(s, "{byte:02x}").unwrap();
    }
    s
}

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod analyze_tests {
    use super::*;
    use crate::interpret::Confidence;
    use crate::mock::{MockAlignmentSource, MockReferenceSource};
    use crate::reader::{AlignmentHeader, SequenceInfo};

    const TEST_MD5: [u8; 16] = [0x12; 16];
    const WRONG_CHR17_LENGTH: u64 = 12_345_678; // not GRCh37 or GRCh38

    fn small_test_config() -> AnalysisConfig {
        AnalysisConfig {
            pmp22_region: Region::new("chr17", 1_000_000, 1_050_000),
            chr2_control_region: Region::new("chr2", 5_000_000, 5_100_000),
            boundary_scan_region: Region::new("chr17", 10_000_000, 10_300_000),
            // Test-sized run length window so the 100 kb elevated sub-region
            // in the CN=3 happy-path test classifies as classical.
            min_classical_run_length: 50_000,
            max_classical_run_length: 200_000,
            ..AnalysisConfig::default()
        }
    }

    fn valid_header(md5: [u8; 16]) -> AlignmentHeader {
        AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: GRCH38_CHR17_LENGTH,
                md5: Some(md5),
            }],
        }
    }

    fn header_with_wrong_length() -> AlignmentHeader {
        AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: WRONG_CHR17_LENGTH,
                md5: Some(TEST_MD5),
            }],
        }
    }

    fn header_missing_chr17() -> AlignmentHeader {
        AlignmentHeader {
            sequences: vec![
                SequenceInfo {
                    name: "chr1".to_string(),
                    length: 248_956_422,
                    md5: None,
                },
                SequenceInfo {
                    name: "chr2".to_string(),
                    length: 242_193_529,
                    md5: None,
                },
            ],
        }
    }

    fn header_with_md5_none() -> AlignmentHeader {
        AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: GRCH38_CHR17_LENGTH,
                md5: None,
            }],
        }
    }

    fn build_mock_cn2(config: &AnalysisConfig, header: AlignmentHeader) -> MockAlignmentSource {
        MockAlignmentSource::with_coverage(&config.pmp22_region, 30)
            .add_coverage(&config.chr2_control_region, 30)
            .add_coverage(&config.boundary_scan_region, 30)
            .with_header(header)
    }

    fn build_mock_cn3_cmt1a(
        config: &AnalysisConfig,
        header: AlignmentHeader,
    ) -> MockAlignmentSource {
        // 100 kb elevated sub-region inside the 300 kb boundary region.
        let elevated = Region::new("chr17", 10_100_000, 10_200_000);
        MockAlignmentSource::with_coverage(&config.pmp22_region, 45)
            .add_coverage(&config.chr2_control_region, 30)
            .add_coverage(&config.boundary_scan_region, 30)
            .add_coverage(&elevated, 15)
            .with_header(header)
    }

    #[test]
    fn analyze_happy_path_cn2() {
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, valid_header(TEST_MD5));
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let result = analyze(&mut source, &mut reference, &config).unwrap();
        assert_eq!(result.depth.estimated_cn, 2);
        assert_eq!(result.interpretation.confidence, Confidence::Full);
        assert_eq!(result.interpretation.subtype_call, None);
        assert!(result.interpretation.atypical_flags.is_empty());
    }

    #[derive(Debug, PartialEq, Eq)]
    enum ProgressKind {
        HeaderVerified,
        WalkStart(&'static str),
        WalkDone(&'static str),
        DepthRatioDone,
        BoundaryScanDone,
        InterpretDone,
    }

    /// Collapse a `Progress` event to its phase tag, or `None` for
    /// `WalkTick` which is sampled rather than a phase boundary.
    fn kind_of(event: &Progress) -> Option<ProgressKind> {
        match event {
            Progress::HeaderVerified => Some(ProgressKind::HeaderVerified),
            Progress::WalkStart { phase, .. } => Some(ProgressKind::WalkStart(phase)),
            Progress::WalkDone { phase, .. } => Some(ProgressKind::WalkDone(phase)),
            Progress::DepthRatioDone { .. } => Some(ProgressKind::DepthRatioDone),
            Progress::BoundaryScanDone { .. } => Some(ProgressKind::BoundaryScanDone),
            Progress::InterpretDone { .. } => Some(ProgressKind::InterpretDone),
            Progress::WalkTick { .. } => None,
        }
    }

    #[test]
    fn analyze_with_progress_emits_expected_event_sequence() {
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, valid_header(TEST_MD5));
        let mut reference = MockReferenceSource::new(TEST_MD5);

        let mut events: Vec<ProgressKind> = Vec::new();
        analyze_with_progress(&mut source, &mut reference, &config, |event| {
            if let Some(kind) = kind_of(&event) {
                events.push(kind);
            }
        })
        .unwrap();

        let expected = vec![
            ProgressKind::HeaderVerified,
            ProgressKind::WalkStart("depth_ratio:pmp22"),
            ProgressKind::WalkDone("depth_ratio:pmp22"),
            ProgressKind::WalkStart("depth_ratio:chr2_control"),
            ProgressKind::WalkDone("depth_ratio:chr2_control"),
            ProgressKind::WalkStart("depth_ratio:rai1"),
            ProgressKind::WalkDone("depth_ratio:rai1"),
            ProgressKind::DepthRatioDone,
            ProgressKind::WalkStart("boundary_scan:chr17"),
            ProgressKind::WalkDone("boundary_scan:chr17"),
            ProgressKind::BoundaryScanDone,
            ProgressKind::InterpretDone,
        ];
        assert_eq!(events, expected);
    }

    #[test]
    fn analyze_with_progress_noop_callback_matches_analyze() {
        // The two entry points must produce byte-identical AnalysisResults
        // on the same inputs.
        let config = small_test_config();
        let mut src1 = build_mock_cn2(&config, valid_header(TEST_MD5));
        let mut ref1 = MockReferenceSource::new(TEST_MD5);
        let direct = analyze(&mut src1, &mut ref1, &config).unwrap();

        let mut src2 = build_mock_cn2(&config, valid_header(TEST_MD5));
        let mut ref2 = MockReferenceSource::new(TEST_MD5);
        let via_progress = analyze_with_progress(&mut src2, &mut ref2, &config, |_| {}).unwrap();

        assert_eq!(direct, via_progress);
    }

    #[test]
    fn analyze_happy_path_cn3_cmt1a() {
        let config = small_test_config();
        let mut source = build_mock_cn3_cmt1a(&config, valid_header(TEST_MD5));
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let result = analyze(&mut source, &mut reference, &config).unwrap();
        assert_eq!(result.depth.estimated_cn, 3);
        assert_eq!(result.interpretation.confidence, Confidence::Full);
        assert_eq!(
            result.interpretation.subtype_call,
            Some(interpret::Subtype::CMT1A)
        );
        assert!(result.interpretation.atypical_flags.is_empty());
    }

    #[test]
    fn analyze_fails_on_md5_mismatch() {
        let config = small_test_config();
        let header_md5 = [0x11; 16];
        let ref_md5 = [0x22; 16];
        let mut source = build_mock_cn2(&config, valid_header(header_md5));
        let mut reference = MockReferenceSource::new(ref_md5);
        let err = match analyze(&mut source, &mut reference, &config) {
            Ok(r) => panic!("expected MD5 mismatch, got {r:?}"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::ReferenceMd5Mismatch { .. }),
            "expected ReferenceMd5Mismatch, got {err:?}"
        );
        let msg = err.to_string();
        assert!(msg.contains("11111111"), "expected hex 11... in: {msg}");
        assert!(msg.contains("22222222"), "expected hex 22... in: {msg}");
        assert_eq!(err.exit_code(), 13);
    }

    #[test]
    fn analyze_fails_on_wrong_build_chr17_length() {
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, header_with_wrong_length());
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let err = match analyze(&mut source, &mut reference, &config) {
            Ok(r) => panic!("expected wrong build, got {r:?}"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::ReferenceBuildMismatch { .. }),
            "expected ReferenceBuildMismatch, got {err:?}"
        );
        let msg = err.to_string();
        assert!(msg.contains("81195210"), "expected GRCh37 length in: {msg}");
        assert!(msg.contains("83257441"), "expected GRCh38 length in: {msg}");
        assert_eq!(err.exit_code(), 12);
    }

    #[test]
    fn analyze_fails_on_missing_chr17_header() {
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, header_missing_chr17());
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let err = match analyze(&mut source, &mut reference, &config) {
            Ok(r) => panic!("expected wrong build, got {r:?}"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::ReferenceBuildMismatch { .. }),
            "expected ReferenceBuildMismatch, got {err:?}"
        );
        assert!(err.to_string().contains("no chr17 or 17 sequence"));
    }

    #[test]
    fn analyze_skips_md5_when_header_has_no_m5_tag() {
        // Header has chr17 at the right length but no M5 tag.
        // The MD5 check is soft-skipped; analyze proceeds normally.
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, header_with_md5_none());
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let result = analyze(&mut source, &mut reference, &config).unwrap();
        assert_eq!(result.depth.estimated_cn, 2);
        assert_eq!(result.interpretation.confidence, Confidence::Full);
    }

    #[test]
    fn analyze_refuses_below_coverage_floor() {
        // Autosomal coverage 10× is below the 15× floor. Interpret returns
        // Refused; analyze propagates that as Ok, not Err.
        let config = small_test_config();
        let mut source = MockAlignmentSource::with_coverage(&config.pmp22_region, 10)
            .add_coverage(&config.chr2_control_region, 10)
            .add_coverage(&config.boundary_scan_region, 10)
            .with_header(valid_header(TEST_MD5));
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let result = analyze(&mut source, &mut reference, &config).unwrap();
        assert_eq!(result.interpretation.confidence, Confidence::Refused);
        assert_eq!(result.interpretation.copy_number, None);
    }

    #[test]
    fn analyze_fails_on_empty_pmp22_region() {
        // Mock has records only in chr2 control and boundary regions;
        // nothing in PMP22. depth_ratio returns EmptyRegion, propagated by
        // analyze.
        let config = small_test_config();
        let mut source = MockAlignmentSource::new(valid_header(TEST_MD5), Vec::new())
            .add_coverage(&config.chr2_control_region, 30)
            .add_coverage(&config.boundary_scan_region, 30);
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let err = match analyze(&mut source, &mut reference, &config) {
            Ok(r) => panic!("expected EmptyRegion, got {r:?}"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::EmptyRegion { .. }),
            "expected EmptyRegion, got {err:?}"
        );
    }

    #[test]
    fn analyze_fails_on_empty_header() {
        // Empty header (zero sequences) should produce a distinct error
        // message pointing at the actual cause - "empty header" - rather
        // than the generic "no chr17 or 17" message.
        let config = small_test_config();
        let empty_header = AlignmentHeader {
            sequences: Vec::new(),
        };
        let mut source = MockAlignmentSource::new(empty_header, Vec::new());
        let mut reference = MockReferenceSource::new(TEST_MD5);
        let err = match analyze(&mut source, &mut reference, &config) {
            Ok(r) => panic!("expected empty-header error, got {r:?}"),
            Err(e) => e,
        };
        assert!(
            matches!(err, Error::ReferenceBuildMismatch { .. }),
            "expected ReferenceBuildMismatch, got {err:?}"
        );
        assert!(
            err.to_string().contains("empty header"),
            "expected 'empty header' in the error message; got: {err}"
        );
    }

    #[test]
    fn analyze_walks_each_region_exactly_once() {
        // analyze calls query_region four times: PMP22, chr2 control, and
        // RAI1 (all from depth_ratio), then the boundary-scan region.
        let config = small_test_config();
        let mut source = build_mock_cn2(&config, valid_header(TEST_MD5));
        let mut reference = MockReferenceSource::new(TEST_MD5);
        analyze(&mut source, &mut reference, &config).unwrap();
        assert_eq!(source.query_count(), 4);
    }

    #[test]
    fn default_config_validates() {
        AnalysisConfig::default().validate().unwrap();
    }

    #[test]
    fn validate_rejects_negative_coverage_floor() {
        let config = AnalysisConfig {
            coverage_floor_refuse: -1.0,
            ..AnalysisConfig::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn validate_rejects_inverted_coverage_floors() {
        let config = AnalysisConfig {
            coverage_floor_low: 10.0, // less than refuse (15.0)
            ..AnalysisConfig::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn validate_rejects_bad_duplication_threshold() {
        let config = AnalysisConfig {
            duplication_threshold: 0.9, // must be > 1.0
            ..AnalysisConfig::default()
        };
        assert!(config.validate().is_err());
    }

    #[test]
    fn validate_rejects_inverted_run_lengths() {
        let config = AnalysisConfig {
            min_classical_run_length: 3_000_000,
            max_classical_run_length: 1_000_000,
            ..AnalysisConfig::default()
        };
        assert!(config.validate().is_err());
    }
}
