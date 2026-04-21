//! Classification decision tree.
//!
//! Combines `DepthRatio` and `BoundaryScanResult` into a final
//! [`Interpretation`]. Pure function - no I/O, no `Result` wrapping. The
//! decision tree is total: every (CN, boundary, coverage) input maps to a
//! defined `Interpretation`.
//!
//! See the interpretation table, the 1-2 Mb classical run range rule, and
//! the 15x / 25x coverage floor for the rules.

use crate::boundary::{BoundaryDirection, BoundaryRun, BoundaryScanResult};
use crate::depth::DepthRatio;
use crate::AnalysisConfig;

/// Confidence tier from the coverage-floor rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Confidence {
    /// ≥ 25× autosomal coverage. Standard report.
    Full,
    /// 15× ≤ coverage < 25×. CN call returned but flagged for independent
    /// confirmation; `plain_language` is prepended with "Low-confidence".
    Low,
    /// < 15× coverage. No CN call returned; `copy_number` is `None`.
    Refused,
}

/// Subtypes this tool can classify in v0.1. Matches the
/// `SubtypeStatus::Supported` rows in [`crate::matrix::SUPPORTED`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Subtype {
    CMT1A,
    HNPP,
    /// Potocki-Lupski syndrome - RAI1 duplication (17p11.2).
    PTLS,
    /// Smith-Magenis syndrome - RAI1 deletion (17p11.2).
    SMS,
    /// Yuan-Harel-Lupski syndrome - contiguous PMP22+RAI1 duplication.
    YUHAL,
}

/// Reasons a finding is marked atypical or unusual. v0.2+ may add variants
/// as new detection methods land.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtypicalFlag {
    /// CN=3 with a duplication run shorter than `min_classical_run_length`.
    Cmt1aBoundaryTooShort,
    /// CN=3 with a duplication run longer than `max_classical_run_length`
    /// (Yuan-Harel-Lupski territory).
    Cmt1aBoundaryTooLong,
    /// CN=3 but no contiguous boundary run found.
    Cmt1aNoBoundaryRun,
    /// CN=1 with a deletion run length outside the classical HNPP range.
    HnppBoundaryAtypical,
    /// CN=1 but no contiguous boundary run found.
    HnppNoBoundaryRun,
    /// Homozygous deletion (CN=0) - very rare, severe HNPP.
    HomozygousDeletion,
    /// Triplication or higher (CN ≥ 4).
    Triplication,
    /// Boundary run direction contradicts the copy-number call (e.g. CN=3
    /// with a deletion run). Strongly suggests low-quality data or a mosaic
    /// finding.
    BoundaryDirectionMismatch,
}

/// Final interpretation returned by [`interpret`].
#[derive(Debug, Clone, PartialEq)]
pub struct Interpretation {
    /// Copy number estimate. `None` when `confidence == Refused`.
    pub copy_number: Option<u8>,
    pub confidence: Confidence,
    /// Called syndrome, or `None` for normal (CN=2), refused (coverage
    /// floor), or atypical cases.
    pub subtype_call: Option<Subtype>,
    /// Short user-visible string ready to embed above the disclaimer in the
    /// text report and the web result panel.
    pub plain_language: String,
    /// Diagnostic flags indicating why a call is atypical. Empty when the
    /// result is a clean CN=2 or a clean CMT1A / HNPP classification.
    pub atypical_flags: Vec<AtypicalFlag>,
}

/// Combine the depth-ratio and boundary-scan results into a final
/// [`Interpretation`].
pub fn interpret(
    depth: &DepthRatio,
    boundary: &BoundaryScanResult,
    config: &AnalysisConfig,
    build: crate::ReferenceBuild,
) -> Interpretation {
    let autosomal = depth.autosomal_mean_depth;

    // Defensive NaN / infinity guard. `NaN < coverage_floor` is false, so
    // without this guard NaN would slip through the floor check and yield a
    // bogus classification.
    if !autosomal.is_finite() {
        return refused_for_non_finite_coverage();
    }

    // Coverage-floor gate.
    if autosomal < config.coverage_floor_refuse {
        return refused(autosomal, config.coverage_floor_refuse);
    }
    let confidence = if autosomal < config.coverage_floor_low {
        Confidence::Low
    } else {
        Confidence::Full
    };

    // Per-CN classification.
    let cn = depth.estimated_cn;
    let (subtype_call, atypical_flags) = match cn {
        0 => (None, vec![AtypicalFlag::HomozygousDeletion]),
        1 => classify_cn1(boundary, config, build),
        2 => (None, vec![]),
        3 => classify_cn3(boundary, config, build),
        _ => (None, vec![AtypicalFlag::Triplication]),
    };

    let plain_language =
        build_plain_language(cn, subtype_call, &atypical_flags, confidence, depth, config);

    Interpretation {
        copy_number: Some(cn),
        confidence,
        subtype_call,
        plain_language,
        atypical_flags,
    }
}

/// Classify a CN=1 (deletion) finding using REP-aware boundary detection,
/// with fallback to length-based classification for test configs.
fn classify_cn1(
    boundary: &BoundaryScanResult,
    config: &AnalysisConfig,
    build: crate::ReferenceBuild,
) -> (Option<Subtype>, Vec<AtypicalFlag>) {
    match boundary.longest_run() {
        None => (None, vec![AtypicalFlag::HnppNoBoundaryRun]),
        Some(run) if run.direction == BoundaryDirection::Duplication => {
            (None, vec![AtypicalFlag::BoundaryDirectionMismatch])
        }
        Some(run) => {
            let (sub, flags) = classify_run_by_rep_region(run, true, build, config);
            if sub.is_some() {
                return (sub, flags);
            }
            if run.length >= config.min_classical_run_length
                && run.length <= config.max_classical_run_length
            {
                (Some(Subtype::HNPP), vec![])
            } else {
                (None, flags)
            }
        }
    }
}

fn classify_cn3(
    boundary: &BoundaryScanResult,
    config: &AnalysisConfig,
    build: crate::ReferenceBuild,
) -> (Option<Subtype>, Vec<AtypicalFlag>) {
    match boundary.longest_run() {
        None => (None, vec![AtypicalFlag::Cmt1aNoBoundaryRun]),
        Some(run) if run.direction == BoundaryDirection::Deletion => {
            (None, vec![AtypicalFlag::BoundaryDirectionMismatch])
        }
        Some(run) => {
            let (sub, flags) = classify_run_by_rep_region(run, false, build, config);
            if sub.is_some() {
                return (sub, flags);
            }
            if run.length >= config.min_classical_run_length
                && run.length <= config.max_classical_run_length
            {
                (Some(Subtype::CMT1A), vec![])
            } else {
                (None, flags)
            }
        }
    }
}

/// Classify a boundary run by checking which REP cluster(s) it aligns with.
/// The tolerance (`AnalysisConfig::rep_match_tolerance`) accounts for 10 kb
/// window resolution and segdup-edge noise.
fn classify_run_by_rep_region(
    run: &BoundaryRun,
    is_deletion: bool,
    build: crate::ReferenceBuild,
    config: &crate::AnalysisConfig,
) -> (Option<Subtype>, Vec<AtypicalFlag>) {
    let ds = match build {
        crate::ReferenceBuild::GRCh37 => crate::reference_data::grch37(),
        crate::ReferenceBuild::GRCh38 => crate::reference_data::grch38(),
    };
    let tolerance = config.rep_match_tolerance;

    let rep = |needle: &str| {
        ds.reps
            .iter()
            .find(|r| r.name.starts_with(needle))
            .unwrap_or_else(|| panic!("reference_data is missing REP '{needle}'"))
    };
    let cmt1a_distal = rep("CMT1A-REP distal");
    let cmt1a_proximal = rep("CMT1A-REP proximal");
    let sms_distal = rep("SMS-REP distal");
    let sms_proximal = rep("SMS-REP proximal");

    let near = |pos: u32, target: u32| -> bool { pos.abs_diff(target) <= tolerance };
    let run_start = run.start;
    let run_end = run.end;

    let in_cmt1a = near(run_start, cmt1a_distal.region.start.value)
        && near(run_end, cmt1a_proximal.region.end.value);
    let in_sms = near(run_start, sms_distal.region.start.value)
        && near(run_end, sms_proximal.region.end.value);
    let is_yuhal = near(run_start, cmt1a_distal.region.start.value)
        && near(run_end, sms_proximal.region.end.value);

    if is_yuhal && !is_deletion {
        return (Some(Subtype::YUHAL), vec![]);
    }

    if in_cmt1a {
        return if is_deletion {
            (Some(Subtype::HNPP), vec![])
        } else {
            (Some(Subtype::CMT1A), vec![])
        };
    }

    if in_sms {
        return if is_deletion {
            (Some(Subtype::SMS), vec![])
        } else {
            (Some(Subtype::PTLS), vec![])
        };
    }

    // Run doesn't align with any known REP cluster.
    let flag = if is_deletion {
        AtypicalFlag::HnppBoundaryAtypical
    } else if run.length < config.atypical_short_threshold {
        AtypicalFlag::Cmt1aBoundaryTooShort
    } else {
        AtypicalFlag::Cmt1aBoundaryTooLong
    };
    (None, vec![flag])
}

fn refused(autosomal: f64, floor: f64) -> Interpretation {
    let plain_language = format!(
        "Autosomal coverage is {autosomal:.1}× which is below the {floor:.1}× minimum required for reliable detection. No copy-number call is returned. Obtain a higher-coverage sample or proceed directly to clinical MLPA."
    );
    Interpretation {
        copy_number: None,
        confidence: Confidence::Refused,
        subtype_call: None,
        plain_language,
        atypical_flags: vec![],
    }
}

/// Defensive path for NaN/infinite autosomal coverage.
fn refused_for_non_finite_coverage() -> Interpretation {
    Interpretation {
        copy_number: None,
        confidence: Confidence::Refused,
        subtype_call: None,
        plain_language:
            "Autosomal coverage is not a finite number (possibly a data-integrity issue upstream). No copy-number call is returned. Please report the input file as a bug."
                .to_string(),
        atypical_flags: vec![],
    }
}

fn build_plain_language(
    cn: u8,
    subtype: Option<Subtype>,
    flags: &[AtypicalFlag],
    confidence: Confidence,
    depth: &DepthRatio,
    config: &AnalysisConfig,
) -> String {
    let core = core_message(cn, subtype, flags, config);
    if confidence == Confidence::Low {
        format!(
            "Low-confidence call (autosomal coverage {:.1}× is below the {:.1}× full-confidence threshold). Independent confirmation is strongly recommended. {core}",
            depth.autosomal_mean_depth, config.coverage_floor_low
        )
    } else {
        core
    }
}

fn core_message(
    cn: u8,
    subtype: Option<Subtype>,
    flags: &[AtypicalFlag],
    config: &AnalysisConfig,
) -> String {
    // CN=0 - homozygous deletion
    if cn == 0 {
        return "Homozygous PMP22 deletion detected - this is consistent with severe HNPP. This is a very rare finding and requires urgent clinical confirmation by MLPA at a certified diagnostic laboratory.".to_string();
    }

    // CN=2 - normal
    if cn == 2 {
        return "Normal copy number at PMP22 and RAI1. No CMT1A duplication, HNPP deletion, PTLS duplication, or SMS deletion detected.".to_string();
    }

    // CN >= 4 - triplication
    if cn >= 4 {
        return format!(
            "Copy number {cn} detected at PMP22 - consistent with a PMP22 triplication or higher. This is a rare severe form of CMT1A. Requires urgent clinical genetics evaluation."
        );
    }

    // CN=1 branches
    if cn == 1 {
        if subtype == Some(Subtype::HNPP) {
            return "PMP22 deletion detected at the classical HNPP locus. This is consistent with HNPP (hereditary neuropathy with liability to pressure palsies). Requires clinical confirmation by MLPA or array-CGH at a certified diagnostic laboratory.".to_string();
        }
        if subtype == Some(Subtype::SMS) {
            return "RAI1 deletion detected at the Smith-Magenis syndrome locus (17p11.2). This is consistent with SMS. Requires clinical confirmation by MLPA or array-CGH at a certified diagnostic laboratory.".to_string();
        }
        if flags.contains(&AtypicalFlag::HnppBoundaryAtypical) {
            return format!(
                "Copy number 1 detected at PMP22 with an atypical boundary (expected {min}-{max} bases for classical HNPP). Requires clinical evaluation by a neurologist or clinical geneticist.",
                min = config.min_classical_run_length,
                max = config.max_classical_run_length,
            );
        }
        if flags.contains(&AtypicalFlag::BoundaryDirectionMismatch) {
            return "Internal inconsistency: depth ratio suggests CN=1 but the boundary scan detected a duplication. This may indicate low-quality data or a mosaic finding. Requires clinical evaluation.".to_string();
        }
        if flags.contains(&AtypicalFlag::HnppNoBoundaryRun) {
            return "Copy number 1 detected at PMP22 but no contiguous boundary run was found. This may indicate a partial or atypical deletion. Requires clinical evaluation.".to_string();
        }
    }

    // CN=3 branches
    if cn == 3 {
        if subtype == Some(Subtype::CMT1A) {
            return "PMP22 duplication detected at the classical CMT1A locus. This is consistent with CMT1A, the most common form of inherited peripheral neuropathy. Requires clinical confirmation by MLPA or array-CGH at a certified diagnostic laboratory.".to_string();
        }
        if subtype == Some(Subtype::PTLS) {
            return "RAI1 duplication detected at the Potocki-Lupski syndrome locus (17p11.2). This is consistent with PTLS. Requires clinical confirmation by MLPA or array-CGH at a certified diagnostic laboratory.".to_string();
        }
        if subtype == Some(Subtype::YUHAL) {
            return "Contiguous PMP22-RAI1 duplication detected spanning the CMT1A and PTLS loci. This is consistent with Yuan-Harel-Lupski syndrome (YUHAL), a rare contiguous gene duplication syndrome. Requires urgent clinical genetics evaluation.".to_string();
        }
        if flags.contains(&AtypicalFlag::Cmt1aBoundaryTooShort) {
            return format!(
                "Copy number 3 detected at PMP22 with a duplication run shorter than the classical CMT1A size (minimum {min} bases). This may indicate a partial duplication, mosaic CMT1A, or a marginal signal. Requires clinical evaluation.",
                min = config.min_classical_run_length,
            );
        }
        if flags.contains(&AtypicalFlag::Cmt1aBoundaryTooLong) {
            return format!(
                "Copy number 3 detected at PMP22 with a duplication run longer than the classical CMT1A size (maximum {max} bases). This may indicate a Yuan-Harel-Lupski-type contiguous gene duplication. Requires clinical genetics evaluation.",
                max = config.max_classical_run_length,
            );
        }
        if flags.contains(&AtypicalFlag::BoundaryDirectionMismatch) {
            return "Internal inconsistency: depth ratio suggests CN=3 but the boundary scan detected a deletion. This may indicate low-quality data or a mosaic finding. Requires clinical evaluation.".to_string();
        }
        if flags.contains(&AtypicalFlag::Cmt1aNoBoundaryRun) {
            return "Copy number 3 detected at PMP22 but no contiguous boundary run was found. This may indicate mosaic CMT1A or a marginal signal. Requires clinical evaluation.".to_string();
        }
    }

    // Soft fallback: a new `AtypicalFlag` variant without a matching branch
    // should degrade to a generic message rather than crash the report.
    format!("Copy number {cn} detected at PMP22. Requires clinical evaluation.")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::{BoundaryDirection, BoundaryRun, BoundaryScanResult, WindowDepth};
    use crate::depth::DepthRatio;

    fn fake_depth(cn: u8, autosomal: f64) -> DepthRatio {
        let pmp22 = autosomal * (cn as f64 / 2.0);
        let ratio = if autosomal > 0.0 {
            pmp22 / autosomal
        } else {
            0.0
        };
        DepthRatio {
            pmp22_mean_depth: pmp22,
            autosomal_mean_depth: autosomal,
            ratio,
            estimated_cn: cn,
            rai1_mean_depth: 0.0,
            rai1_ratio: 0.0,
            rai1_estimated_cn: 0,
        }
    }

    fn fake_boundary(run: Option<BoundaryRun>) -> BoundaryScanResult {
        BoundaryScanResult {
            windows: Vec::<WindowDepth>::new(),
            runs: run.into_iter().collect(),
        }
    }

    fn classical_dup_run() -> BoundaryRun {
        BoundaryRun {
            start: 14_190_000,
            end: 15_580_000,
            length: 1_390_001,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        }
    }

    fn classical_del_run() -> BoundaryRun {
        BoundaryRun {
            start: 14_190_000,
            end: 15_580_000,
            length: 1_390_001,
            mean_normalized_depth: 0.5,
            direction: BoundaryDirection::Deletion,
            bridged_gaps: 0,
        }
    }

    fn dup_run_with_length(length: u32) -> BoundaryRun {
        BoundaryRun {
            start: 14_500_000,
            end: 14_500_000 + length - 1,
            length,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        }
    }

    fn del_run_with_length(length: u32) -> BoundaryRun {
        BoundaryRun {
            start: 14_500_000,
            end: 14_500_000 + length - 1,
            length,
            mean_normalized_depth: 0.5,
            direction: BoundaryDirection::Deletion,
            bridged_gaps: 0,
        }
    }

    // ---------- Coverage band tests ----------

    #[test]
    fn interpret_below_floor_refuses() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 10.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Refused);
        assert_eq!(result.copy_number, None);
        assert_eq!(result.subtype_call, None);
        assert!(result.atypical_flags.is_empty());
        assert!(result.plain_language.contains("10.0×"));
        assert!(result.plain_language.contains("below the 15.0× minimum"));
    }

    #[test]
    fn interpret_at_exactly_15x_is_low() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 15.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Low);
        assert_eq!(result.copy_number, Some(2));
    }

    #[test]
    fn interpret_at_exactly_25x_is_full() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 25.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Full);
        assert_eq!(result.copy_number, Some(2));
    }

    #[test]
    fn interpret_low_confidence_string_prepended() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 20.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Low);
        assert!(result.plain_language.starts_with("Low-confidence call"));
        assert!(result.plain_language.contains("20.0×"));
        // Still includes the standard message at the end.
        assert!(result.plain_language.contains("Normal copy number"));
    }

    // ---------- CN=0 ----------

    #[test]
    fn interpret_cn0_flags_homozygous_deletion() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(0, 30.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Full);
        assert_eq!(result.copy_number, Some(0));
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::HomozygousDeletion]
        );
        assert!(result.plain_language.contains("Homozygous PMP22 deletion"));
    }

    // ---------- CN=1 ----------

    #[test]
    fn interpret_cn1_with_classical_deletion_calls_hnpp() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(1, 30.0);
        let boundary = fake_boundary(Some(classical_del_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, Some(Subtype::HNPP));
        assert!(result.atypical_flags.is_empty());
        assert!(result.plain_language.contains("HNPP"));
    }

    #[test]
    fn interpret_cn1_with_short_deletion_flags_atypical() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(1, 30.0);
        let boundary = fake_boundary(Some(del_run_with_length(500_000))); // < 1 Mb
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::HnppBoundaryAtypical]
        );
    }

    #[test]
    fn interpret_cn1_with_duplication_flags_direction_mismatch() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(1, 30.0);
        let boundary = fake_boundary(Some(classical_dup_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::BoundaryDirectionMismatch]
        );
        assert!(result.plain_language.contains("Internal inconsistency"));
    }

    #[test]
    fn interpret_cn1_with_no_run_flags_no_boundary_run() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(1, 30.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(result.atypical_flags, vec![AtypicalFlag::HnppNoBoundaryRun]);
    }

    // ---------- CN=2 ----------

    #[test]
    fn interpret_cn2_clean_returns_no_flags() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 30.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert!(result.atypical_flags.is_empty());
        assert_eq!(
            result.plain_language,
            "Normal copy number at PMP22 and RAI1. No CMT1A duplication, HNPP deletion, PTLS duplication, or SMS deletion detected."
        );
    }

    #[test]
    fn interpret_cn2_with_spurious_run_still_clean() {
        // CN=2 with a short duplication run (sub-classical noise) must NOT
        // be flagged — CN=2 is always clean.
        let config = AnalysisConfig::default();
        let depth = fake_depth(2, 30.0);
        let boundary = fake_boundary(Some(dup_run_with_length(20_000)));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert!(
            result.atypical_flags.is_empty(),
            "CN=2 should never produce atypical flags, got {:?}",
            result.atypical_flags
        );
    }

    // ---------- CN=3 ----------

    #[test]
    fn interpret_cn3_with_classical_duplication_calls_cmt1a() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(classical_dup_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, Some(Subtype::CMT1A));
        assert!(result.atypical_flags.is_empty());
        assert!(result.plain_language.contains("CMT1A"));
    }

    #[test]
    fn interpret_cn3_with_unaligned_short_dup_flags_atypical() {
        // Run at 14.5M-14.8M - doesn't align with any known REP cluster.
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(dup_run_with_length(300_000))); // < 500 kb
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::Cmt1aBoundaryTooShort]
        );
    }

    #[test]
    fn interpret_cn3_with_unaligned_long_dup_flags_atypical() {
        // Run at 14.5M-17.5M - doesn't align cleanly with any single REP pair.
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(dup_run_with_length(3_000_000))); // > 500 kb, no REP match
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::Cmt1aBoundaryTooLong]
        );
    }

    #[test]
    fn interpret_cn3_with_deletion_flags_direction_mismatch() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(classical_del_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::BoundaryDirectionMismatch]
        );
    }

    #[test]
    fn interpret_cn3_with_no_run_flags_no_boundary_run() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, None);
        assert_eq!(
            result.atypical_flags,
            vec![AtypicalFlag::Cmt1aNoBoundaryRun]
        );
    }

    // ---------- CN ≥ 4 ----------

    #[test]
    fn interpret_refuses_on_non_finite_autosomal() {
        let config = AnalysisConfig::default();
        let depth = DepthRatio {
            pmp22_mean_depth: 30.0,
            autosomal_mean_depth: f64::NAN,
            ratio: f64::NAN,
            estimated_cn: 2,
            rai1_mean_depth: 0.0,
            rai1_ratio: 0.0,
            rai1_estimated_cn: 0,
        };
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.confidence, Confidence::Refused);
        assert_eq!(result.copy_number, None);
        assert_eq!(result.subtype_call, None);
        assert!(result.atypical_flags.is_empty());
        assert!(
            result.plain_language.contains("not a finite number"),
            "expected 'not a finite number' in plain_language, got: {}",
            result.plain_language
        );
    }

    #[test]
    fn interpret_cn4_flags_triplication() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(4, 30.0);
        let boundary = fake_boundary(None);
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.copy_number, Some(4));
        assert_eq!(result.subtype_call, None);
        assert_eq!(result.atypical_flags, vec![AtypicalFlag::Triplication]);
        assert!(result.plain_language.contains("triplication"));
    }

    // --- PTLS / SMS / YUHAL classification ---

    fn ptls_dup_run() -> BoundaryRun {
        // Duplication in the SMS-REP region (GRCh38 UCSC: distal ~16.69M, proximal ~20.56M)
        BoundaryRun {
            start: 16_690_000,
            end: 20_560_000,
            length: 3_870_001,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        }
    }

    fn sms_del_run() -> BoundaryRun {
        BoundaryRun {
            start: 16_690_000,
            end: 20_560_000,
            length: 3_870_001,
            mean_normalized_depth: 0.5,
            direction: BoundaryDirection::Deletion,
            bridged_gaps: 0,
        }
    }

    fn yuhal_dup_run() -> BoundaryRun {
        // Large duplication spanning CMT1A-REP distal through SMS-REP proximal
        BoundaryRun {
            start: 14_170_000,
            end: 20_560_000,
            length: 6_390_001,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        }
    }

    fn atypical_dup_run() -> BoundaryRun {
        // Duplication that doesn't align with any known REP cluster
        BoundaryRun {
            start: 16_000_000,
            end: 16_400_000,
            length: 400_001,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        }
    }

    #[test]
    fn classify_ptls_duplication() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(ptls_dup_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, Some(Subtype::PTLS), "expected PTLS");
    }

    #[test]
    fn classify_sms_deletion() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(1, 30.0);
        let boundary = fake_boundary(Some(sms_del_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, Some(Subtype::SMS), "expected SMS");
    }

    #[test]
    fn classify_yuhal_duplication() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(yuhal_dup_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(result.subtype_call, Some(Subtype::YUHAL), "expected YUHAL");
    }

    #[test]
    fn classify_atypical_not_misclassified() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        let boundary = fake_boundary(Some(atypical_dup_run()));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh38);
        assert_eq!(
            result.subtype_call, None,
            "atypical rearrangement should not classify as any known syndrome"
        );
        assert!(
            !result.atypical_flags.is_empty(),
            "should have atypical flags"
        );
    }

    #[test]
    fn classify_cmt1a_with_grch37_build() {
        let config = AnalysisConfig::default();
        let depth = fake_depth(3, 30.0);
        // GRCh37 CMT1A-REP positions: distal ~14.07M, proximal ~15.49M
        let run = BoundaryRun {
            start: 14_074_000,
            end: 15_495_000,
            length: 1_421_001,
            mean_normalized_depth: 1.5,
            direction: BoundaryDirection::Duplication,
            bridged_gaps: 0,
        };
        let boundary = fake_boundary(Some(run));
        let result = interpret(&depth, &boundary, &config, crate::ReferenceBuild::GRCh37);
        assert_eq!(
            result.subtype_call,
            Some(Subtype::CMT1A),
            "GRCh37 CMT1A should classify correctly with GRCh37 build"
        );
    }
}
