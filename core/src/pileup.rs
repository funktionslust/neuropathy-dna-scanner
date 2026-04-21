//! Targeted pileup genotyping for known-pathogenic variants.
//!
//! For each variant in the catalog, queries the alignment at that
//! exact position and counts ref vs alt alleles among high-quality
//! reads. This is the SNV/indel complement to the structural variant
//! detection in `depth.rs` and `boundary.rs`.

use crate::catalog::{ClinSig, Inheritance, VariantCatalog};
use serde::Serialize;

/// Calling thresholds for pileup genotyping.
#[derive(Debug, Clone)]
pub struct PileupConfig {
    /// Minimum number of reads supporting the alt allele. Default: 3.
    pub min_alt_reads: u32,
    /// Minimum allele balance for a heterozygous call. Default: 0.2.
    pub min_het_ab: f64,
    /// Maximum allele balance for a heterozygous call. Default: 0.8.
    pub max_het_ab: f64,
    /// Minimum MAPQ for a read to contribute to pileup. Default: 20.
    pub min_mapq: u8,
    /// Minimum Phred base quality at the variant position. Bases below
    /// this quality are excluded from counting. Default: 20.
    pub min_base_qual: u8,
    /// Minimum total depth at a position. Default: 10.
    pub min_depth: u32,
}

impl Default for PileupConfig {
    fn default() -> Self {
        Self {
            min_alt_reads: 3,
            min_het_ab: 0.2,
            max_het_ab: 0.8,
            min_mapq: 20,
            min_base_qual: 20,
            min_depth: 10,
        }
    }
}

/// Result of checking a single variant position.
#[derive(Debug, Clone, Serialize)]
pub enum VariantCall {
    /// No reads cover this position.
    NoCoverage,
    /// All reads show the reference allele.
    ReferenceOnly { depth: u32 },
    /// Heterozygous: significant alt allele reads present.
    Heterozygous {
        ref_count: u32,
        alt_count: u32,
        allele_balance: f64,
    },
    /// Homozygous alt: nearly all reads show alt.
    Homozygous {
        ref_count: u32,
        alt_count: u32,
        allele_balance: f64,
    },
    /// Some alt reads present but below calling thresholds.
    Marginal {
        ref_count: u32,
        alt_count: u32,
        allele_balance: f64,
    },
    /// Alt reads present but all on one strand - likely artifact.
    StrandBias {
        ref_count: u32,
        alt_count: u32,
        allele_balance: f64,
        forward_alt: u32,
        reverse_alt: u32,
    },
    /// Insufficient depth at this position.
    LowCoverage { depth: u32 },
    /// Hemizygous: present on a single-copy chromosome (chrX/Y in males).
    Hemizygous { ref_count: u32, alt_count: u32 },
}

impl VariantCall {
    /// Whether this call found the variant (het, hom, or hemizygous).
    pub fn is_positive(&self) -> bool {
        matches!(
            self,
            VariantCall::Heterozygous { .. }
                | VariantCall::Homozygous { .. }
                | VariantCall::Hemizygous { .. }
        )
    }
}

/// A checked variant: the catalog entry plus the pileup result.
#[derive(Debug, Clone, Serialize)]
pub struct CheckedVariant {
    pub gene: String,
    pub chrom: String,
    pub pos: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub clinvar_id: u32,
    pub significance: ClinSig,
    pub review_stars: u8,
    pub condition: String,
    pub inheritance: Inheritance,
    pub consequence: String,
    pub call: VariantCall,
}

/// Summary of pileup screening results.
#[derive(Debug, Clone, Serialize)]
pub struct PileupResult {
    /// Variants where the alt allele was detected (het or hom).
    pub positive_findings: Vec<CheckedVariant>,
    /// Total variants checked.
    pub total_checked: usize,
    /// Variants with no coverage.
    pub no_coverage_count: usize,
    /// Variants with low coverage (below min_depth).
    pub low_coverage_count: usize,
    /// Catalog metadata (date, gene count, etc.)
    pub catalog_date: String,
    pub catalog_gene_count: usize,
    pub catalog_variant_count: usize,
}

/// Strand-aware variant classification. Rejects calls where all alt
/// reads are on a single strand (forward or reverse only), which is
/// a classic library-prep / alignment artifact.
pub fn classify_pileup_with_strand(
    ref_count: u32,
    alt_count: u32,
    forward_alt: u32,
    reverse_alt: u32,
    config: &PileupConfig,
) -> VariantCall {
    let base_call = classify_pileup(ref_count, alt_count, config);
    match &base_call {
        VariantCall::Heterozygous { .. } | VariantCall::Homozygous { .. } => {
            if forward_alt == 0 || reverse_alt == 0 {
                let total = ref_count + alt_count;
                VariantCall::StrandBias {
                    ref_count,
                    alt_count,
                    allele_balance: if total > 0 {
                        alt_count as f64 / total as f64
                    } else {
                        0.0
                    },
                    forward_alt,
                    reverse_alt,
                }
            } else {
                base_call
            }
        }
        _ => base_call,
    }
}

/// Sex-aware variant classification. Adjusts min_depth for hemizygous
/// chromosomes (chrX/Y in males) and calls hemizygous genotypes.
pub fn classify_pileup_with_context(
    ref_count: u32,
    alt_count: u32,
    forward_alt: u32,
    reverse_alt: u32,
    config: &PileupConfig,
    chrom: &str,
    sex: crate::Sex,
) -> VariantCall {
    let is_hemi = matches!(
        (chrom, sex),
        ("chrX" | "X", crate::Sex::Male) | ("chrY" | "Y", crate::Sex::Male)
    );

    if is_hemi {
        // Hemizygous: halved depth threshold, haploid genotype calling
        let total = ref_count + alt_count;
        let effective_min = config.min_depth / 2;
        if total == 0 {
            return VariantCall::NoCoverage;
        }
        if total < effective_min {
            return VariantCall::LowCoverage { depth: total };
        }
        if alt_count == 0 {
            return VariantCall::ReferenceOnly { depth: total };
        }
        let ab = alt_count as f64 / total as f64;
        if ab >= 0.8 {
            VariantCall::Hemizygous {
                ref_count,
                alt_count,
            }
        } else if ab < 0.2 {
            VariantCall::ReferenceOnly { depth: total }
        } else {
            VariantCall::Marginal {
                ref_count,
                alt_count,
                allele_balance: ab,
            }
        }
    } else {
        // Diploid: use strand-aware classification
        classify_pileup_with_strand(ref_count, alt_count, forward_alt, reverse_alt, config)
    }
}

/// Classify pileup counts into a variant call.
pub fn classify_pileup(ref_count: u32, alt_count: u32, config: &PileupConfig) -> VariantCall {
    let total = ref_count + alt_count;

    if total == 0 {
        return VariantCall::NoCoverage;
    }

    if total < config.min_depth {
        return VariantCall::LowCoverage { depth: total };
    }

    if alt_count == 0 {
        return VariantCall::ReferenceOnly { depth: total };
    }

    let ab = alt_count as f64 / total as f64;

    if alt_count >= config.min_alt_reads && ab >= config.max_het_ab {
        VariantCall::Homozygous {
            ref_count,
            alt_count,
            allele_balance: ab,
        }
    } else if alt_count >= config.min_alt_reads
        && ab >= config.min_het_ab
        && ab <= config.max_het_ab
    {
        VariantCall::Heterozygous {
            ref_count,
            alt_count,
            allele_balance: ab,
        }
    } else {
        VariantCall::Marginal {
            ref_count,
            alt_count,
            allele_balance: ab,
        }
    }
}

/// Build a PileupResult from a list of checked variants.
pub fn build_pileup_result(checked: Vec<CheckedVariant>, catalog: &VariantCatalog) -> PileupResult {
    let total_checked = checked.len();
    let no_coverage_count = checked
        .iter()
        .filter(|c| matches!(c.call, VariantCall::NoCoverage))
        .count();
    let low_coverage_count = checked
        .iter()
        .filter(|c| matches!(c.call, VariantCall::LowCoverage { .. }))
        .count();
    let positive_findings: Vec<CheckedVariant> = checked
        .into_iter()
        .filter(|c| c.call.is_positive())
        .collect();

    PileupResult {
        positive_findings,
        total_checked,
        no_coverage_count,
        low_coverage_count,
        catalog_date: catalog.meta.generated_date.clone(),
        catalog_gene_count: catalog.meta.gene_count,
        catalog_variant_count: catalog.meta.variant_count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_config() -> PileupConfig {
        PileupConfig::default()
    }

    #[test]
    fn no_coverage_returns_no_coverage() {
        let call = classify_pileup(0, 0, &default_config());
        assert!(matches!(call, VariantCall::NoCoverage));
    }

    #[test]
    fn ref_only_returns_ref_only() {
        let call = classify_pileup(30, 0, &default_config());
        assert!(matches!(call, VariantCall::ReferenceOnly { depth: 30 }));
    }

    #[test]
    fn het_call_with_good_balance() {
        let call = classify_pileup(15, 12, &default_config());
        assert!(matches!(call, VariantCall::Heterozygous { .. }));
        assert!(call.is_positive());
    }

    #[test]
    fn hom_call_with_high_alt() {
        let call = classify_pileup(2, 28, &default_config());
        assert!(matches!(call, VariantCall::Homozygous { .. }));
        assert!(call.is_positive());
    }

    #[test]
    fn marginal_with_few_alt_reads() {
        // 2 alt reads, below min_alt_reads=3
        let call = classify_pileup(28, 2, &default_config());
        assert!(matches!(call, VariantCall::Marginal { .. }));
        assert!(!call.is_positive());
    }

    #[test]
    fn low_coverage_below_min_depth() {
        let call = classify_pileup(4, 3, &default_config());
        assert!(matches!(call, VariantCall::LowCoverage { depth: 7 }));
    }

    #[test]
    fn marginal_with_low_allele_balance() {
        // 3 alt reads but AB = 3/30 = 0.1, below min_het_ab=0.2
        let call = classify_pileup(27, 3, &default_config());
        assert!(matches!(call, VariantCall::Marginal { .. }));
    }
}
