//! Plain-text report renderer.
//!
//! Produces the report the CLI writes next to the input file and the
//! web app offers for download. Contract: `render_report` is pure - given
//! identical inputs it produces byte-identical output. The only
//! source of non-determinism is the `timestamp` field in `ReportMeta`, which
//! the caller provides.

use std::fmt::Write as _;

use crate::boundary::BoundaryDirection;
use crate::matrix::{self, DISCLAIMER};
use crate::AnalysisResult;

const SECTION_RULE_DOUBLE: &str =
    "================================================================";
const SECTION_RULE_SINGLE: &str =
    "----------------------------------------------------------------";
const WRAP_WIDTH: usize = 72;

/// Metadata about the run, injected by the caller (CLI driver or WASM
/// frontend). Kept separate from [`AnalysisResult`] so the core algorithm
/// layer never has to know about filenames or timestamps.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportMeta {
    /// Tool version string, typically `env!("CARGO_PKG_VERSION")`.
    pub tool_version: &'static str,
    /// ISO 8601 / RFC 3339 timestamp in UTC (e.g. `"2026-04-15T14:23:47Z"`).
    /// Caller is responsible for formatting.
    pub timestamp: String,
    /// Input file basename for display. Full path is intentionally omitted
    /// to avoid leaking user directory structure in reports.
    pub input_filename: String,
    /// chr17 MD5 hex string from the CRAM header, if present.
    pub input_md5_hex: Option<String>,
}

/// Render the full analysis result as a plain-text report.
pub fn render_report(result: &AnalysisResult, meta: &ReportMeta) -> String {
    let mut out = String::with_capacity(4096);

    // Section 1: matrix (single source of truth).
    out.push_str(&matrix::render_text());
    out.push('\n');

    // Section 2: metadata.
    writeln!(out, "{SECTION_RULE_DOUBLE}").unwrap();
    writeln!(
        out,
        "neuropathy-dna-scanner v{} -- analysis report",
        meta.tool_version
    )
    .unwrap();
    writeln!(out, "{SECTION_RULE_DOUBLE}").unwrap();
    out.push('\n');
    writeln!(out, "Timestamp:       {}", meta.timestamp).unwrap();
    writeln!(out, "Input:           {}", meta.input_filename).unwrap();
    match &meta.input_md5_hex {
        Some(md5) => writeln!(out, "Reference MD5:   {md5}").unwrap(),
        None => writeln!(out, "Reference MD5:   not present in header").unwrap(),
    }
    out.push('\n');

    // Section 3: depth-ratio analysis.
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    writeln!(out, "Depth-ratio analysis").unwrap();
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    writeln!(
        out,
        "PMP22 region mean depth:  {:.1}×",
        result.depth.pmp22_mean_depth
    )
    .unwrap();
    writeln!(
        out,
        "Autosomal mean depth:     {:.1}×",
        result.depth.autosomal_mean_depth
    )
    .unwrap();
    writeln!(out, "Ratio:                    {:.3}", result.depth.ratio).unwrap();
    writeln!(
        out,
        "Estimated copy number:    {}",
        result.depth.estimated_cn
    )
    .unwrap();
    out.push('\n');

    // Section 4: boundary scan.
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    writeln!(out, "Boundary scan").unwrap();
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    match result.boundary.longest_run() {
        Some(run) => {
            let dir = match run.direction {
                BoundaryDirection::Duplication => "Duplication",
                BoundaryDirection::Deletion => "Deletion",
            };
            writeln!(
                out,
                "Longest run:           chr17:{}-{}",
                run.start, run.end
            )
            .unwrap();
            writeln!(out, "Direction:             {dir}").unwrap();
            writeln!(out, "Length:                {} bases", run.length).unwrap();
            writeln!(
                out,
                "Mean normalized depth: {:.3}",
                run.mean_normalized_depth
            )
            .unwrap();
        }
        None => {
            writeln!(out, "No contiguous duplicated or deleted run detected.").unwrap();
        }
    }
    out.push('\n');

    // Section 5: interpretation.
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    writeln!(out, "Interpretation").unwrap();
    writeln!(out, "{SECTION_RULE_SINGLE}").unwrap();
    match result.interpretation.copy_number {
        Some(cn) => writeln!(out, "Copy number:    {cn}").unwrap(),
        None => writeln!(out, "Copy number:    (refused - coverage below floor)").unwrap(),
    }
    writeln!(
        out,
        "Confidence:     {:?}",
        result.interpretation.confidence
    )
    .unwrap();
    match result.interpretation.subtype_call {
        Some(subtype) => writeln!(out, "Subtype:        {subtype:?}").unwrap(),
        None => writeln!(out, "Subtype:        none").unwrap(),
    }
    let flags_str = if result.interpretation.atypical_flags.is_empty() {
        "none".to_string()
    } else {
        result
            .interpretation
            .atypical_flags
            .iter()
            .map(|f| format!("{f:?}"))
            .collect::<Vec<_>>()
            .join(", ")
    };
    writeln!(out, "Atypical flags: {flags_str}").unwrap();
    out.push('\n');
    out.push_str(&wrap_at_72(&result.interpretation.plain_language));
    out.push_str("\n\n");

    // Section 6: disclaimer.
    writeln!(out, "{SECTION_RULE_DOUBLE}").unwrap();
    out.push_str(DISCLAIMER);
    out.push('\n');
    writeln!(out, "{SECTION_RULE_DOUBLE}").unwrap();

    out
}

/// Break `text` into lines of at most [`WRAP_WIDTH`] characters on word
/// boundaries. Preserves existing `\n` characters (splits on them first,
/// wraps each paragraph, rejoins). Single words longer than `WRAP_WIDTH`
/// are emitted unwrapped.
fn wrap_at_72(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + 16);
    let mut first_paragraph = true;
    for paragraph in text.split('\n') {
        if !first_paragraph {
            out.push('\n');
        }
        first_paragraph = false;
        let mut line_len = 0usize;
        let mut first_word = true;
        for word in paragraph.split_whitespace() {
            let word_len = word.chars().count();
            if first_word {
                out.push_str(word);
                line_len = word_len;
                first_word = false;
            } else if line_len + 1 + word_len <= WRAP_WIDTH {
                out.push(' ');
                out.push_str(word);
                line_len += 1 + word_len;
            } else {
                out.push('\n');
                out.push_str(word);
                line_len = word_len;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::{BoundaryDirection, BoundaryRun, BoundaryScanResult};
    use crate::depth::DepthRatio;
    use crate::interpret::{Confidence, Interpretation, Subtype};

    fn fake_meta() -> ReportMeta {
        ReportMeta {
            tool_version: "0.1.0",
            timestamp: "2026-04-15T14:23:47Z".to_string(),
            input_filename: "sample.cram".to_string(),
            input_md5_hex: Some("1d3085cf3e275ecb01ed2b30db8b6e2c".to_string()),
        }
    }

    fn build_result_cn2() -> AnalysisResult {
        AnalysisResult {
            depth: DepthRatio {
                pmp22_mean_depth: 30.0,
                autosomal_mean_depth: 30.0,
                ratio: 1.0,
                estimated_cn: 2, rai1_mean_depth: 0.0, rai1_ratio: 0.0, rai1_estimated_cn: 0,
            },
            boundary: BoundaryScanResult {
                windows: Vec::new(),
                runs: vec![],
            },
            interpretation: Interpretation {
                copy_number: Some(2),
                confidence: Confidence::Full,
                subtype_call: None,
                plain_language:
                    "Normal copy number at PMP22 and RAI1. No CMT1A duplication, HNPP deletion, PTLS duplication, or SMS deletion detected."
                        .to_string(),
                atypical_flags: vec![],
            },
        }
    }

    fn build_result_cn3_cmt1a() -> AnalysisResult {
        AnalysisResult {
            depth: DepthRatio {
                pmp22_mean_depth: 45.0,
                autosomal_mean_depth: 30.0,
                ratio: 1.5,
                estimated_cn: 3, rai1_mean_depth: 0.0, rai1_ratio: 0.0, rai1_estimated_cn: 0,
            },
            boundary: BoundaryScanResult {
                windows: Vec::new(),
                runs: vec![BoundaryRun {
                    start: 14_190_000,
                    end: 15_580_000,
                    length: 1_390_001,
                    mean_normalized_depth: 1.5,
                    direction: BoundaryDirection::Duplication,
                    bridged_gaps: 0,
                }],
            },
            interpretation: Interpretation {
                copy_number: Some(3),
                confidence: Confidence::Full,
                subtype_call: Some(Subtype::CMT1A),
                plain_language:
                    "PMP22 duplication detected at the classical CMT1A locus. This is consistent with CMT1A, the most common form of inherited peripheral neuropathy. Requires clinical confirmation by MLPA or array-CGH at a certified diagnostic laboratory."
                        .to_string(),
                atypical_flags: vec![],
            },
        }
    }

    fn build_result_refused() -> AnalysisResult {
        AnalysisResult {
            depth: DepthRatio {
                pmp22_mean_depth: 10.0,
                autosomal_mean_depth: 10.0,
                ratio: 1.0,
                estimated_cn: 2, rai1_mean_depth: 0.0, rai1_ratio: 0.0, rai1_estimated_cn: 0,
            },
            boundary: BoundaryScanResult {
                windows: Vec::new(),
                runs: vec![],
            },
            interpretation: Interpretation {
                copy_number: None,
                confidence: Confidence::Refused,
                subtype_call: None,
                plain_language:
                    "Autosomal coverage is 10.0× which is below the 15.0× minimum required for reliable detection. No copy-number call is returned. Obtain a higher-coverage sample or proceed directly to clinical MLPA."
                        .to_string(),
                atypical_flags: vec![],
            },
        }
    }

    // ---------- wrap_at_72 tests ----------

    #[test]
    fn wrap_at_72_short_text_is_unchanged() {
        let input = "short line under 72 chars";
        assert_eq!(wrap_at_72(input), "short line under 72 chars");
    }

    #[test]
    fn wrap_at_72_long_text_breaks_on_word_boundary() {
        let input = "This is a much longer sentence that absolutely exceeds the seventy-two character limit and needs to wrap.";
        let wrapped = wrap_at_72(input);
        for line in wrapped.split('\n') {
            assert!(
                line.chars().count() <= WRAP_WIDTH,
                "line {line:?} exceeds {WRAP_WIDTH}"
            );
        }
        // The wrapped output should contain at least one newline.
        assert!(wrapped.contains('\n'));
    }

    #[test]
    fn wrap_at_72_preserves_explicit_newlines() {
        let input = "First paragraph.\nSecond paragraph.";
        let wrapped = wrap_at_72(input);
        assert_eq!(wrapped, "First paragraph.\nSecond paragraph.");
    }

    #[test]
    fn wrap_at_72_empty_string_returns_empty() {
        assert_eq!(wrap_at_72(""), "");
    }

    // Snapshot tests. Regenerate by running with NDS_UPDATE_SNAPSHOTS=1 set
    // in the environment.

    const EXPECTED_CN2: &str = include_str!("../tests/snapshots/report_cn2.txt");
    const EXPECTED_CN3_CMT1A: &str = include_str!("../tests/snapshots/report_cn3_cmt1a.txt");
    const EXPECTED_REFUSED: &str = include_str!("../tests/snapshots/report_refused.txt");

    fn assert_or_update_snapshot(rel_path: &str, actual: &str, expected: &str) {
        if std::env::var("NDS_UPDATE_SNAPSHOTS").is_ok() {
            let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(rel_path);
            std::fs::write(path, actual).expect("write snapshot");
            return;
        }
        assert_eq!(actual, expected);
    }

    #[test]
    fn snapshot_cn2_normal_report() {
        let report = render_report(&build_result_cn2(), &fake_meta());
        assert_or_update_snapshot("tests/snapshots/report_cn2.txt", &report, EXPECTED_CN2);
    }

    #[test]
    fn snapshot_cn3_cmt1a_report() {
        let report = render_report(&build_result_cn3_cmt1a(), &fake_meta());
        assert_or_update_snapshot(
            "tests/snapshots/report_cn3_cmt1a.txt",
            &report,
            EXPECTED_CN3_CMT1A,
        );
    }

    #[test]
    fn snapshot_refused_report() {
        let report = render_report(&build_result_refused(), &fake_meta());
        assert_or_update_snapshot(
            "tests/snapshots/report_refused.txt",
            &report,
            EXPECTED_REFUSED,
        );
    }

    #[test]
    fn render_report_is_deterministic() {
        let result = build_result_cn2();
        let meta = fake_meta();
        let first = render_report(&result, &meta);
        let second = render_report(&result, &meta);
        assert_eq!(first, second);
    }

    #[test]
    fn render_report_handles_missing_md5_header() {
        // Exercises the `input_md5_hex: None` branch, which the full
        // snapshot tests don't cover.
        let meta = ReportMeta {
            tool_version: "0.1.0",
            timestamp: "2026-04-15T14:23:47Z".to_string(),
            input_filename: "sample.bam".to_string(),
            input_md5_hex: None,
        };
        let report = render_report(&build_result_cn2(), &meta);
        assert!(report.contains("Reference MD5:   not present in header"));
        assert!(report.contains("Normal copy number at PMP22"));
        assert!(report.contains("Research-grade screening tool only."));
    }
}
