//! Embedded variant catalog for targeted SNV/indel screening.
//!
//! The catalog is built offline from ClinVar + PanelApp data by the
//! `build_catalog` tool and embedded in the binary at compile time.
//! Every variant is traceable to its ClinVar accession number.

use serde::{Deserialize, Serialize};

/// Clinical significance from ClinVar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ClinSig {
    Pathogenic,
    LikelyPathogenic,
}

/// Mode of inheritance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Inheritance {
    AutosomalDominant,
    AutosomalRecessive,
    XLinked,
    /// Both AD and AR reported for this gene.
    Both,
    Unknown,
}

/// A single known-pathogenic variant from ClinVar.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathogenicVariant {
    pub chrom: String,
    /// 1-based position.
    pub pos: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub gene: String,
    /// ClinVar variation ID (e.g., 12345).
    pub clinvar_id: u32,
    pub significance: ClinSig,
    /// ClinVar review status stars (0-4).
    pub review_stars: u8,
    pub condition: String,
    pub inheritance: Inheritance,
    /// Molecular consequence (missense, nonsense, frameshift, splice, etc.).
    pub consequence: String,
}

/// Metadata for the embedded catalog.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CatalogMeta {
    /// Date the catalog was generated.
    pub generated_date: String,
    /// ClinVar release date or file date.
    pub clinvar_date: String,
    /// PanelApp panel ID and version.
    pub panelapp_panel: String,
    /// Number of genes in the panel.
    pub gene_count: usize,
    /// Number of variants in the catalog.
    pub variant_count: usize,
    /// Minimum ClinVar review stars included.
    pub min_review_stars: u8,
    /// Genes excluded from the catalog with reasons.
    #[serde(default)]
    pub excluded_genes: Vec<ExcludedGene>,
}

/// A gene excluded from screening with the reason and recommendation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExcludedGene {
    pub symbol: String,
    pub reason: ExclusionReason,
    pub recommendation: String,
}

/// Why a gene was excluded from the catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExclusionReason {
    /// Near-identical paralog prevents reliable short-read genotyping (e.g. SMN1/SMN2).
    ParalogMappability,
    /// Mitochondrial variants require heteroplasmy-aware calling.
    MitochondrialHeteroplasmy,
    /// Segmental duplication causes unreliable depth signals.
    SegmentalDuplication,
    /// Insufficient validation data for this methodology.
    InsufficientValidation,
}

/// The full variant catalog: metadata + variant list.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantCatalog {
    pub meta: CatalogMeta,
    pub variants: Vec<PathogenicVariant>,
}

impl VariantCatalog {
    /// Load a catalog from JSON bytes.
    pub fn from_json(data: &[u8]) -> Result<Self, String> {
        serde_json::from_slice(data).map_err(|e| format!("catalog parse error: {e}"))
    }

    /// Variants for a specific chromosome, sorted by position.
    pub fn variants_on_chrom(&self, chrom: &str) -> Vec<&PathogenicVariant> {
        let mut vs: Vec<&PathogenicVariant> =
            self.variants.iter().filter(|v| v.chrom == chrom).collect();
        vs.sort_by_key(|v| v.pos);
        vs
    }

    /// Distinct gene symbols in the catalog.
    pub fn genes(&self) -> Vec<&str> {
        let mut genes: Vec<&str> = self.variants.iter().map(|v| v.gene.as_str()).collect();
        genes.sort_unstable();
        genes.dedup();
        genes
    }
}
