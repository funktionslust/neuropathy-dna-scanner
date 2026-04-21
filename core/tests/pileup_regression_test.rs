//! Regression tests for pileup classification.
//!
//! These tests exercise the classification functions directly with
//! inputs that match the fixture spike parameters. They guard against
//! regressions in:
//! - Strand bias detection (PNPLA6 false positive)
//! - Sex-aware hemizygous calling (chrX in males)
//! - Multi-alt same-position discrimination

use nds_core::pileup::{
    classify_pileup, classify_pileup_with_context, classify_pileup_with_strand, PileupConfig,
    VariantCall,
};
use nds_core::Sex;

fn default_config() -> PileupConfig {
    PileupConfig::default()
}

// =========================================================================
// Strand bias regression (PNPLA6 false positive)
// =========================================================================

#[test]
fn strand_bias_all_reverse_rejects_het_call() {
    // Simulates PNPLA6 false positive: 7 ref, 3 alt but all alt on reverse
    let call = classify_pileup_with_strand(7, 3, 0, 3, &default_config());
    assert!(
        matches!(call, VariantCall::StrandBias { .. }),
        "expected StrandBias, got {call:?}"
    );
    assert!(!call.is_positive(), "strand bias should not be positive");
}

#[test]
fn strand_bias_all_forward_rejects_het_call() {
    let call = classify_pileup_with_strand(15, 12, 12, 0, &default_config());
    assert!(
        matches!(call, VariantCall::StrandBias { .. }),
        "expected StrandBias, got {call:?}"
    );
}

#[test]
fn strand_balanced_allows_het_call() {
    // Both strands represented — should pass strand check
    let call = classify_pileup_with_strand(15, 12, 6, 6, &default_config());
    assert!(
        matches!(call, VariantCall::Heterozygous { .. }),
        "expected Heterozygous, got {call:?}"
    );
    assert!(call.is_positive());
}

#[test]
fn strand_bias_does_not_affect_ref_only() {
    // No alt reads — strand bias check is irrelevant
    let call = classify_pileup_with_strand(30, 0, 0, 0, &default_config());
    assert!(matches!(call, VariantCall::ReferenceOnly { .. }));
}

#[test]
fn strand_bias_does_not_affect_marginal() {
    // Below min_alt_reads — classify as marginal, not strand bias
    let call = classify_pileup_with_strand(28, 2, 0, 2, &default_config());
    assert!(
        matches!(call, VariantCall::Marginal { .. }),
        "expected Marginal (below min_alt_reads), got {call:?}"
    );
}

// =========================================================================
// Sex-aware classification (chrX hemizygous)
// =========================================================================

#[test]
fn male_chrx_hemizygous_alt_detected() {
    // Male chrX: hemizygous, most reads show alt
    let call = classify_pileup_with_context(1, 14, 7, 7, &default_config(), "chrX", Sex::Male);
    assert!(
        matches!(call, VariantCall::Hemizygous { .. }),
        "expected Hemizygous, got {call:?}"
    );
    assert!(call.is_positive());
}

#[test]
fn male_chrx_reference_only() {
    // Male chrX: all ref, no alt
    let call = classify_pileup_with_context(15, 0, 0, 0, &default_config(), "chrX", Sex::Male);
    assert!(
        matches!(call, VariantCall::ReferenceOnly { .. }),
        "expected ReferenceOnly, got {call:?}"
    );
}

#[test]
fn male_chrx_low_ab_is_reference() {
    // Male chrX: low alt fraction (below 0.2) → treated as ref
    let call = classify_pileup_with_context(14, 1, 1, 0, &default_config(), "chrX", Sex::Male);
    assert!(
        matches!(call, VariantCall::ReferenceOnly { .. }),
        "expected ReferenceOnly for low AB on hemizygous chrX, got {call:?}"
    );
}

#[test]
fn male_chrx_intermediate_ab_is_marginal() {
    // Male chrX: AB between 0.2 and 0.8 on hemizygous — marginal
    let call = classify_pileup_with_context(8, 7, 3, 4, &default_config(), "chrX", Sex::Male);
    assert!(
        matches!(call, VariantCall::Marginal { .. }),
        "expected Marginal for intermediate AB on hemizygous chrX, got {call:?}"
    );
}

#[test]
fn male_chrx_halved_depth_threshold() {
    // Male chrX: depth threshold is halved (min_depth/2 = 5)
    // Total depth 4 < 5 → LowCoverage
    let call = classify_pileup_with_context(3, 1, 1, 0, &default_config(), "chrX", Sex::Male);
    assert!(
        matches!(call, VariantCall::LowCoverage { .. }),
        "expected LowCoverage for depth 4 on hemizygous chrX (threshold=5), got {call:?}"
    );
}

#[test]
fn female_chrx_uses_diploid_calling() {
    // Female chrX: diploid, same as autosomal
    let call = classify_pileup_with_context(15, 12, 6, 6, &default_config(), "chrX", Sex::Female);
    assert!(
        matches!(call, VariantCall::Heterozygous { .. }),
        "expected Heterozygous (diploid) for female chrX, got {call:?}"
    );
}

#[test]
fn autosomal_ignores_sex() {
    // Autosomal chromosome: sex doesn't matter, uses standard diploid calling
    let call = classify_pileup_with_context(15, 12, 6, 6, &default_config(), "chr17", Sex::Male);
    assert!(
        matches!(call, VariantCall::Heterozygous { .. }),
        "expected Heterozygous for autosomal chr17 regardless of sex, got {call:?}"
    );
}

// =========================================================================
// Multi-alt same-position discrimination
// =========================================================================

#[test]
fn correct_alt_detected_at_shared_position() {
    // At a position with C>A and C>T in catalog:
    // If reads show C>A (ref=44, alt=47), classify as het
    let call = classify_pileup(44, 47, &default_config());
    assert!(matches!(call, VariantCall::Heterozygous { .. }));
}

#[test]
fn wrong_alt_not_detected_at_shared_position() {
    // Same position, but checking C>T: ref=44, alt=0 → ref only
    let call = classify_pileup(44, 0, &default_config());
    assert!(matches!(call, VariantCall::ReferenceOnly { .. }));
}

// =========================================================================
// Base quality threshold
// =========================================================================

#[test]
fn config_has_min_base_qual_20() {
    let config = PileupConfig::default();
    assert_eq!(
        config.min_base_qual, 20,
        "min_base_qual default must be 20 to filter PNPLA6-like artifacts"
    );
}

#[test]
fn config_has_min_mapq_20() {
    let config = PileupConfig::default();
    assert_eq!(config.min_mapq, 20);
}
