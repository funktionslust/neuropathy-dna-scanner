//! nds-wasm: WebAssembly bindings for the nds-core library.
//!
//! Production exports:
//! - `version()` - core library version string.
//! - `render_matrix_text()` / `render_matrix_html()` - supported-subtypes matrix.
//! - `required_ranges_bam()` - BAI -> file byte ranges for sparse fetching.
//! - `analyze_bam()` - full analysis pipeline over prefetched BAM + BAI bytes.

mod pileup_wasm;

use serde::Serialize;
use wasm_bindgen::prelude::*;

use nds_core::catalog::VariantCatalog;
use nds_core::pileup::PileupConfig;
use nds_core::reader::{AlignmentHeader, Region};
use nds_core::reference::ReferenceSource;
use nds_core::{
    analyze_with_progress, AlignmentSource, AnalysisConfig, CursorBamSource, CursorCramSource,
    Progress, Result as CoreResult,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_sex(s: &str) -> nds_core::Sex {
    match s.to_lowercase().as_str() {
        "male" | "m" => nds_core::Sex::Male,
        "female" | "f" => nds_core::Sex::Female,
        _ => nds_core::Sex::Unknown,
    }
}

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------

#[wasm_bindgen(start)]
pub fn _init() {
    console_error_panic_hook::set_once();
}

// ---------------------------------------------------------------------------
// Simple exports
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn version() -> String {
    nds_core::version().to_string()
}

#[wasm_bindgen]
pub fn render_matrix_text() -> String {
    nds_core::matrix::render_text()
}

#[wasm_bindgen]
pub fn render_matrix_html() -> String {
    nds_core::matrix::render_html()
}

// ---------------------------------------------------------------------------
// Sex detection from CRAI + header
// ---------------------------------------------------------------------------

/// Estimate biological sex from CRAI alignment spans per chromosome.
/// Uses chrX/chrY span ratios relative to autosomal average as a proxy
/// for depth ratios. Returns "male", "female", or "unknown".
#[wasm_bindgen]
pub fn detect_sex_from_index(crai_bytes: &[u8], header_bytes: &[u8]) -> String {
    let result = (|| -> Result<String, String> {
        let index: Vec<noodles_cram::crai::Record> =
            noodles_cram::crai::io::Reader::new(std::io::Cursor::new(crai_bytes))
                .read_index()
                .map_err(|e| format!("CRAI parse: {e}"))?;

        let mut cram_reader = noodles_cram::io::reader::Builder::default()
            .build_from_reader(std::io::Cursor::new(header_bytes));
        let header = cram_reader
            .read_header()
            .map_err(|e| format!("header parse: {e}"))?;

        // Sum slice_length (compressed bytes) per reference sequence as a
        // proxy for read count. This is more reliable than alignment_span
        // which can span the whole chromosome even with few reads.
        let mut bytes_by_ref: std::collections::HashMap<usize, u64> =
            std::collections::HashMap::new();
        for rec in &index {
            if let Some(ref_id) = rec.reference_sequence_id() {
                *bytes_by_ref.entry(ref_id).or_default() += rec.slice_length();
            }
        }

        // Normalize by chromosome length to get a depth proxy.
        let refs = header.reference_sequences();
        let mut chrx_depth: f64 = 0.0;
        let mut chry_depth: f64 = 0.0;
        let mut auto_depths: Vec<f64> = Vec::new();

        for (idx, (name, map)) in refs.iter().enumerate() {
            let name_str = std::str::from_utf8(name.as_ref()).unwrap_or("");
            let bytes = bytes_by_ref.get(&idx).copied().unwrap_or(0);
            let chrom_len = map.length().get() as f64;
            if chrom_len < 1.0 {
                continue;
            }
            let depth_proxy = bytes as f64 / chrom_len;
            match name_str {
                "chrX" | "X" => chrx_depth = depth_proxy,
                "chrY" | "Y" => chry_depth = depth_proxy,
                n if depth_proxy > 0.0
                    && !n.contains('M')
                    && !n.contains("Un")
                    && !n.contains("random")
                    && !n.contains("hap")
                    && !n.contains("alt")
                    && !n.contains("decoy") =>
                {
                    auto_depths.push(depth_proxy);
                }
                _ => {}
            }
        }

        if auto_depths.is_empty() {
            return Ok("unknown".into());
        }

        let auto_mean = auto_depths.iter().sum::<f64>() / auto_depths.len() as f64;
        let sex = nds_core::detect_sex(chrx_depth, chry_depth, auto_mean);
        Ok(match sex {
            nds_core::Sex::Male => "male",
            nds_core::Sex::Female => "female",
            nds_core::Sex::Unknown => "unknown",
        }
        .into())
    })();
    result.unwrap_or_else(|_| "unknown".into())
}

// ---------------------------------------------------------------------------
// Reference data
// ---------------------------------------------------------------------------

/// Returns the full reference dataset (coordinates, genes, REPs, syndromes)
/// as JSON for the given genome build. Defaults to GRCh38.
#[wasm_bindgen]
pub fn reference_data_json(build: &str) -> String {
    let dataset = match build {
        "GRCh37" | "hg19" => nds_core::reference_data::grch37(),
        _ => nds_core::reference_data::grch38(),
    };
    serde_json::to_string_pretty(&dataset).unwrap_or_else(|_| "{}".to_string())
}

// ---------------------------------------------------------------------------
// Variant catalog
// ---------------------------------------------------------------------------

/// Returns the full screening manifest: everything this tool checks for.
/// Combines SV syndromes from reference_data and SNV conditions from the
/// catalog. This is the single source of truth for the UI's checklist.
#[wasm_bindgen]
pub fn screening_manifest(build: &str) -> String {
    let ds = match build {
        "GRCh37" | "hg19" => nds_core::reference_data::grch37(),
        _ => nds_core::reference_data::grch38(),
    };
    let catalog = load_catalog();

    // SV syndromes from reference_data
    let sv_items: Vec<serde_json::Value> = ds
        .syndromes
        .iter()
        .map(|s| {
            // PTLS and SMS are developmental-delay syndromes at the same
            // chr17 locus; they fall out of the same depth scan for free
            // but aren't neuropathies themselves, so the UI flags them as
            // extras. YUHAL is a contiguous PMP22+RAI1 duplication, so it
            // causes a CMT1A-type neuropathy plus PTLS-type developmental
            // features - it stays in core scope because it includes the
            // neuropathy component.
            let is_extra = matches!(s.short_name, "PTLS" | "SMS");
            serde_json::json!({
                "type": "sv",
                "short_name": s.short_name,
                "name": s.name,
                "gene": s.mediated_by,
                "mechanism": s.mechanism,
                "omim": s.omim,
                "omim_url": s.omim_url,
                "icd10": s.icd10,
                "icd11": s.icd11,
                "orphanet": s.orphanet,
                "orphanet_url": s.orphanet_url,
                "genereviews_url": s.genereviews_url,
                "extra_check": is_extra,
            })
        })
        .collect();

    // SNV gene groups from catalog - group variants by gene, deduplicate conditions
    let mut gene_map: std::collections::BTreeMap<&str, Vec<&str>> =
        std::collections::BTreeMap::new();
    let mut gene_variant_count: std::collections::HashMap<&str, usize> =
        std::collections::HashMap::new();
    let mut gene_inheritance: std::collections::HashMap<&str, &str> =
        std::collections::HashMap::new();
    for v in &catalog.variants {
        gene_map
            .entry(v.gene.as_str())
            .or_default()
            .push(v.condition.as_str());
        *gene_variant_count.entry(v.gene.as_str()).or_default() += 1;
        gene_inheritance
            .entry(v.gene.as_str())
            .or_insert(match v.inheritance {
                nds_core::catalog::Inheritance::AutosomalDominant => "AD",
                nds_core::catalog::Inheritance::AutosomalRecessive => "AR",
                nds_core::catalog::Inheritance::XLinked => "XL",
                nds_core::catalog::Inheritance::Both => "AD/AR",
                nds_core::catalog::Inheritance::Unknown => "unknown",
            });
    }

    let build_gene_item = |gene: &&str, conditions: &Vec<&str>| {
        // ClinVar conditions are pipe-separated (e.g. "CMT type 2A|not provided").
        // Split on pipes, deduplicate, filter junk, pick the best one.
        let mut all_conds: Vec<String> = conditions
            .iter()
            .flat_map(|c| c.split('|'))
            .map(|s| s.trim().replace('_', " "))
            .filter(|s| !s.is_empty() && s != "not provided" && s != "not specified")
            .collect();
        all_conds.sort_unstable();
        all_conds.dedup();
        // Prefer conditions mentioning the gene name or "Charcot"
        let primary_condition = all_conds
            .iter()
            .find(|c| {
                let l = c.to_lowercase();
                l.contains("charcot")
                    || l.contains(&gene.to_lowercase())
                    || l.contains("neuropathy")
            })
            .or_else(|| all_conds.first())
            .map(|s| s.as_str())
            .unwrap_or("Hereditary neuropathy");
        let is_cmt = all_conds.iter().any(|c| {
            let lower = c.to_lowercase();
            lower.contains("charcot")
                || lower.contains("cmt")
                || lower.contains("hereditary motor and sensory neuropathy")
        });
        // Get primary chromosome for this gene
        let gene_chrom = catalog
            .variants
            .iter()
            .find(|v| v.gene.as_str() == *gene)
            .map(|v| v.chrom.as_str())
            .unwrap_or("");
        serde_json::json!({
            "gene": gene,
            "chrom": gene_chrom,
            "condition": primary_condition,
            "variant_count": gene_variant_count.get(gene).unwrap_or(&0),
            "inheritance": gene_inheritance.get(gene).unwrap_or(&"unknown"),
            "is_cmt_core": is_cmt,
        })
    };

    let mut cmt_core: Vec<serde_json::Value> = Vec::new();
    let mut extended_panel: Vec<serde_json::Value> = Vec::new();
    for (gene, conditions) in &gene_map {
        let item = build_gene_item(gene, conditions);
        if item["is_cmt_core"].as_bool().unwrap_or(false) {
            cmt_core.push(item);
        } else {
            extended_panel.push(item);
        }
    }

    // Count SNV-only variants (single-base ref+alt) vs total
    let snv_count = catalog
        .variants
        .iter()
        .filter(|v| v.ref_allele.len() == 1 && v.alt_allele.len() == 1)
        .count();

    let manifest = serde_json::json!({
        "sv_screening": sv_items,
        "cmt_core_genes": cmt_core,
        "extended_panel_genes": extended_panel,
        "catalog_date": catalog.meta.generated_date,
        "catalog_variant_count": catalog.meta.variant_count,
        "catalog_snv_count": snv_count,
        "catalog_skipped_indels": catalog.meta.variant_count - snv_count,
        "catalog_gene_count": catalog.meta.gene_count,
        "catalog_filter": "Pathogenic + Likely Pathogenic, review status >= 1 star",
        "catalog_source": "ClinVar VCF (GRCh38) filtered by PanelApp panel 85 (Hereditary neuropathy)",
    });

    serde_json::to_string(&manifest).unwrap_or_else(|_| "{}".to_string())
}

/// Returns the full variant catalog as CSV for download/audit.
/// Each row is one ClinVar variant with its ID, position, gene, etc.
#[wasm_bindgen]
pub fn catalog_csv() -> String {
    let catalog = load_catalog();
    let mut csv = String::from(
        "ClinVar_ID,Gene,Chromosome,Position,Ref,Alt,Significance,Review_Stars,Condition,Inheritance,Consequence,ClinVar_URL\n"
    );
    for v in &catalog.variants {
        let sig = match v.significance {
            nds_core::catalog::ClinSig::Pathogenic => "Pathogenic",
            nds_core::catalog::ClinSig::LikelyPathogenic => "Likely_pathogenic",
        };
        let inh = match v.inheritance {
            nds_core::catalog::Inheritance::AutosomalDominant => "AD",
            nds_core::catalog::Inheritance::AutosomalRecessive => "AR",
            nds_core::catalog::Inheritance::XLinked => "XL",
            nds_core::catalog::Inheritance::Both => "AD/AR",
            nds_core::catalog::Inheritance::Unknown => "unknown",
        };
        // Clean condition: take first pipe-separated value
        let condition = v.condition.split('|').next().unwrap_or("").replace('"', "");
        csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},\"{}\",{},{},https://www.ncbi.nlm.nih.gov/clinvar/variation/{}/\n",
            v.clinvar_id, v.gene, v.chrom, v.pos, v.ref_allele, v.alt_allele,
            sig, v.review_stars, condition, inh, v.consequence, v.clinvar_id
        ));
    }
    csv
}

/// Returns the embedded variant catalog metadata as JSON.
#[wasm_bindgen]
pub fn catalog_info() -> String {
    let catalog = load_catalog();
    serde_json::to_string(&catalog.meta).unwrap_or_else(|_| "{}".to_string())
}

/// Returns the number of chromosomes with catalog variants, for range planning.
#[wasm_bindgen]
pub fn catalog_chromosomes() -> String {
    let catalog = load_catalog();
    let mut chrom_set: Vec<&str> = catalog.variants.iter().map(|v| v.chrom.as_str()).collect();
    chrom_set.sort_unstable();
    chrom_set.dedup();
    serde_json::to_string(&chrom_set).unwrap_or_else(|_| "[]".to_string())
}

/// Returns compact genomic regions that cover all catalog variant positions.
/// Each region is a gene-level cluster with a small padding. Used by the
/// browser to compute File.slice() byte ranges for the pileup step.
/// Returns JSON: `[{"chrom":"1","start":11989000,"end":11990000}, ...]`
#[wasm_bindgen]
pub fn catalog_regions() -> String {
    let catalog = load_catalog();
    let regions = catalog_to_regions(&catalog);
    serde_json::to_string(&regions).unwrap_or_else(|_| "[]".to_string())
}

/// Compute compact regions from catalog variant positions.
/// Groups variants by chromosome, then clusters nearby positions
/// (within 10 kb of each other) into single regions with 1 kb padding.
fn catalog_to_regions(catalog: &VariantCatalog) -> Vec<serde_json::Value> {
    use std::collections::BTreeMap;

    // Group positions by chromosome
    let mut by_chrom: BTreeMap<&str, Vec<u32>> = BTreeMap::new();
    for v in &catalog.variants {
        // Only SNVs (single-base ref and alt) for now
        if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
            by_chrom.entry(v.chrom.as_str()).or_default().push(v.pos);
        }
    }

    let mut regions = Vec::new();
    let padding: u32 = 1_000;
    let merge_distance: u32 = 10_000;

    for (chrom, mut positions) in by_chrom {
        positions.sort_unstable();
        positions.dedup();
        if positions.is_empty() {
            continue;
        }

        // Cluster nearby positions into regions
        let mut cluster_start = positions[0].saturating_sub(padding);
        let mut cluster_end = positions[0] + padding;

        for &pos in &positions[1..] {
            if pos <= cluster_end + merge_distance {
                // Extend current cluster
                cluster_end = pos + padding;
            } else {
                // Emit current cluster, start new one
                regions.push(serde_json::json!({
                    "chrom": chrom,
                    "start": cluster_start,
                    "end": cluster_end,
                }));
                cluster_start = pos.saturating_sub(padding);
                cluster_end = pos + padding;
            }
        }
        // Emit last cluster
        regions.push(serde_json::json!({
            "chrom": chrom,
            "start": cluster_start,
            "end": cluster_end,
        }));
    }

    regions
}

/// Run pileup genotyping on a BAM file against the embedded catalog.
/// Returns JSON with positive findings, total checked, coverage stats.
#[wasm_bindgen]
pub fn check_variants_bam(bam_bytes: &[u8], bai_bytes: &[u8]) -> JsValue {
    let catalog = load_catalog();
    let config = PileupConfig::default();

    match pileup_wasm::pileup_bam(bam_bytes, bai_bytes, &catalog, &config) {
        Ok(result) => serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL),
        Err(e) => {
            let error = serde_json::json!({ "status": "error", "error": e });
            serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL)
        }
    }
}

/// Run pileup genotyping on a CRAM file against the embedded catalog.
/// Requires reference chromosome sequences for CRAM decoding.
#[wasm_bindgen]
pub fn check_variants_cram(
    cram_bytes: &[u8],
    crai_bytes: &[u8],
    ref_chr17: &[u8],
    ref_chr2: &[u8],
    sex_str: &str,
) -> JsValue {
    let catalog = load_catalog();
    let config = PileupConfig::default();
    let sex = parse_sex(sex_str);

    // Build reference sequence entries for all chromosomes the catalog needs.
    // For now we have chr17 and chr2; variants on other chromosomes will
    // show as NoCoverage (the CRAM won't have those chroms loaded anyway).
    let ref_seqs = vec![
        ("chr17".to_string(), ref_chr17.to_vec()),
        ("chr2".to_string(), ref_chr2.to_vec()),
    ];

    match pileup_wasm::pileup_cram(cram_bytes, crai_bytes, ref_seqs, &catalog, &config, sex) {
        Ok(result) => serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL),
        Err(e) => {
            let error = serde_json::json!({ "status": "error", "error": e });
            serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL)
        }
    }
}

/// Sparse BAM pileup from File.slice() data.
#[wasm_bindgen]
pub fn check_variants_bam_sparse(data: Vec<u8>, ranges_json: &str, bai_bytes: Vec<u8>) -> JsValue {
    let catalog = load_catalog();
    let config = PileupConfig::default();

    let range_descs: Vec<SparseRange> = match serde_json::from_str(ranges_json) {
        Ok(r) => r,
        Err(e) => {
            let error =
                serde_json::json!({ "status": "error", "error": format!("invalid ranges: {e}") });
            return serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL);
        }
    };
    let descs: Vec<(u64, u64)> = range_descs
        .iter()
        .map(|r| (r.file_offset, r.length))
        .collect();

    match pileup_wasm::pileup_bam_sparse(data, &descs, bai_bytes, &catalog, &config) {
        Ok(result) => serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL),
        Err(e) => {
            let error = serde_json::json!({ "status": "error", "error": e });
            serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL)
        }
    }
}

/// Per-chromosome CRAM pileup. The browser calls this in a loop for
/// each chromosome, loading reference + CRAM data one at a time.
#[wasm_bindgen]
pub fn check_variants_cram_chrom(
    cram_data: Vec<u8>,
    ranges_json: &str,
    crai_bytes: Vec<u8>,
    ref_seq: Vec<u8>,
    chrom_name: &str,
    sex_str: &str,
) -> JsValue {
    let catalog = load_catalog();
    let config = PileupConfig::default();
    let sex = parse_sex(sex_str);

    let range_descs: Vec<SparseRange> = match serde_json::from_str(ranges_json) {
        Ok(r) => r,
        Err(e) => {
            let error =
                serde_json::json!({ "status": "error", "error": format!("invalid ranges: {e}") });
            return serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL);
        }
    };
    let descs: Vec<(u64, u64)> = range_descs
        .iter()
        .map(|r| (r.file_offset, r.length))
        .collect();

    let ref_seqs = vec![(chrom_name.to_string(), ref_seq)];

    // Filter catalog to just this chromosome
    let chrom_catalog = filter_catalog_to_chrom(&catalog, chrom_name);

    match pileup_wasm::pileup_cram_sparse(
        cram_data,
        &descs,
        crai_bytes,
        ref_seqs,
        &chrom_catalog,
        &config,
        sex,
    ) {
        Ok(result) => serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL),
        Err(e) => {
            let error = serde_json::json!({ "status": "error", "error": e });
            serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL)
        }
    }
}

/// Filter the catalog to variants on a single chromosome.
/// Handles both bare ("1") and prefixed ("chr1") naming.
fn filter_catalog_to_chrom(catalog: &VariantCatalog, chrom: &str) -> VariantCatalog {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let with_chr = if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{chrom}")
    };

    let variants: Vec<_> = catalog
        .variants
        .iter()
        .filter(|v| v.chrom == bare || v.chrom == with_chr)
        .cloned()
        .collect();

    VariantCatalog {
        meta: catalog.meta.clone(),
        variants,
    }
}

/// Sparse CRAM pileup from File.slice() data.
#[wasm_bindgen]
pub fn check_variants_cram_sparse(
    data: Vec<u8>,
    ranges_json: &str,
    crai_bytes: Vec<u8>,
    ref_chr17: Vec<u8>,
    ref_chr2: Vec<u8>,
    sex_str: &str,
) -> JsValue {
    let catalog = load_catalog();
    let config = PileupConfig::default();
    let sex = parse_sex(sex_str);

    let range_descs: Vec<SparseRange> = match serde_json::from_str(ranges_json) {
        Ok(r) => r,
        Err(e) => {
            let error =
                serde_json::json!({ "status": "error", "error": format!("invalid ranges: {e}") });
            return serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL);
        }
    };
    let descs: Vec<(u64, u64)> = range_descs
        .iter()
        .map(|r| (r.file_offset, r.length))
        .collect();

    let ref_seqs = vec![
        ("chr17".to_string(), ref_chr17),
        ("chr2".to_string(), ref_chr2),
    ];

    match pileup_wasm::pileup_cram_sparse(
        data, &descs, crai_bytes, ref_seqs, &catalog, &config, sex,
    ) {
        Ok(result) => serde_wasm_bindgen::to_value(&result).unwrap_or(JsValue::NULL),
        Err(e) => {
            let error = serde_json::json!({ "status": "error", "error": e });
            serde_wasm_bindgen::to_value(&error).unwrap_or(JsValue::NULL)
        }
    }
}

fn load_catalog() -> VariantCatalog {
    // `catalog.json.gz` is the git-tracked source of truth; `build.rs`
    // decompresses it into `OUT_DIR/catalog.json` so a fresh checkout
    // builds without a manual pre-step.
    static CATALOG_JSON: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/catalog.json"));
    VariantCatalog::from_json(CATALOG_JSON).expect("embedded catalog is invalid")
}

// ---------------------------------------------------------------------------
// Build detection
// ---------------------------------------------------------------------------

/// Detect the reference build (GRCh37 or GRCh38) from a BAM/CRAM header.
/// Pass the first ~4 MB of the file. Returns "GRCh38", "GRCh37", or an error.
#[wasm_bindgen]
pub fn detect_build(header_bytes: &[u8], is_cram: bool) -> String {
    let result = if is_cram {
        let mut reader = noodles_cram::io::reader::Builder::default()
            .build_from_reader(std::io::Cursor::new(header_bytes));
        reader
            .read_header()
            .map(|h| nds_core::sources::sam_header_to_alignment_header(&h))
    } else {
        let mut reader = noodles_bam::io::Reader::new(std::io::Cursor::new(header_bytes));
        reader
            .read_header()
            .map(|h| nds_core::sources::sam_header_to_alignment_header(&h))
    };

    match result {
        Ok(header) => {
            let chr17 = header.sequences.iter().find(|s| {
                s.name.eq_ignore_ascii_case("chr17") || s.name.eq_ignore_ascii_case("17")
            });
            match chr17 {
                Some(seq) if seq.length == 83_257_441 => "GRCh38".to_string(),
                Some(seq) if seq.length == 81_195_210 => "GRCh37".to_string(),
                Some(seq) => format!("unknown (chr17 length: {})", seq.length),
                None => "unknown (chr17 not found)".to_string(),
            }
        }
        Err(e) => format!("error: {e}"),
    }
}

// ---------------------------------------------------------------------------
// Range query
// ---------------------------------------------------------------------------

#[wasm_bindgen]
pub fn required_ranges_cram(crai_bytes: &[u8], header_bytes: &[u8], regions_json: &str) -> JsValue {
    let reply = match compute_ranges_cram(crai_bytes, header_bytes, regions_json, None) {
        Ok(ranges) => RangesReply::Ok { ranges },
        Err(error) => RangesReply::Error { error },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

/// Tight-margin CRAM range query for SNV pileup (±50 kb instead of ±5 Mb).
#[wasm_bindgen]
pub fn required_ranges_cram_tight(
    crai_bytes: &[u8],
    header_bytes: &[u8],
    regions_json: &str,
) -> JsValue {
    let margin = nds_core::sources::range_query::CRAM_MARGIN_TIGHT;
    let reply = match compute_ranges_cram(crai_bytes, header_bytes, regions_json, Some(margin)) {
        Ok(ranges) => RangesReply::Ok { ranges },
        Err(error) => RangesReply::Error { error },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

fn compute_ranges_cram(
    crai_bytes: &[u8],
    header_bytes: &[u8],
    regions_json: &str,
    margin: Option<u32>,
) -> Result<Vec<RangeDto>, String> {
    let inputs: Vec<RegionInput> =
        serde_json::from_str(regions_json).map_err(|e| format!("invalid regions JSON: {e}"))?;
    let regions: Vec<Region> = inputs
        .iter()
        .map(|r| Region::new(&r.chrom, r.start, r.end))
        .collect();
    let ranges = match margin {
        Some(m) => nds_core::sources::range_query::required_ranges_cram_with_margin(
            crai_bytes,
            header_bytes,
            &regions,
            m,
        ),
        None => {
            nds_core::sources::range_query::required_ranges_cram(crai_bytes, header_bytes, &regions)
        }
    }
    .map_err(|e| format!("CRAM range query failed: {e}"))?;
    Ok(ranges
        .iter()
        .map(|r| RangeDto {
            offset: r.offset,
            length: r.length,
        })
        .collect())
}

#[wasm_bindgen]
pub fn required_ranges_bam(bai_bytes: &[u8], header_bytes: &[u8], regions_json: &str) -> JsValue {
    let reply = match compute_ranges(bai_bytes, header_bytes, regions_json) {
        Ok(ranges) => RangesReply::Ok { ranges },
        Err(error) => RangesReply::Error { error },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

#[derive(Serialize)]
#[serde(tag = "status")]
enum RangesReply {
    #[serde(rename = "ok")]
    Ok { ranges: Vec<RangeDto> },
    #[serde(rename = "error")]
    Error { error: String },
}

fn compute_ranges(
    bai_bytes: &[u8],
    header_bytes: &[u8],
    regions_json: &str,
) -> Result<Vec<RangeDto>, String> {
    let inputs: Vec<RegionInput> =
        serde_json::from_str(regions_json).map_err(|e| format!("invalid regions JSON: {e}"))?;
    let regions: Vec<Region> = inputs
        .iter()
        .map(|r| Region::new(&r.chrom, r.start, r.end))
        .collect();
    let ranges =
        nds_core::sources::range_query::required_ranges_bam(bai_bytes, header_bytes, &regions)
            .map_err(|e| format!("range query failed: {e}"))?;
    Ok(ranges
        .iter()
        .map(|r| RangeDto {
            offset: r.offset,
            length: r.length,
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Full analysis
// ---------------------------------------------------------------------------

/// Run the full analysis pipeline on prefetched BAM + BAI bytes.
///
/// `config_json`: optional JSON with overridable fields:
/// `{"chr2_control_start", "chr2_control_end", "boundary_scan_start", "boundary_scan_end"}`.
/// Missing fields use `AnalysisConfig::default()`.
#[wasm_bindgen]
pub fn analyze_bam(bam_bytes: Vec<u8>, bai_bytes: Vec<u8>, config_json: &str) -> JsValue {
    let start = js_sys::Date::now() as u64;
    let reply = match run_analyze(bam_bytes, bai_bytes, config_json) {
        Ok(mut ok) => {
            ok.elapsed_ms = (js_sys::Date::now() as u64).saturating_sub(start);
            AnalyzeReply::Ok(ok)
        }
        Err(message) => AnalyzeReply::Error { message },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

/// Run the full analysis pipeline on prefetched CRAM + CRAI bytes.
/// No reference FASTA needed - only alignment positions are read.
#[wasm_bindgen]
pub fn analyze_cram(cram_bytes: Vec<u8>, crai_bytes: Vec<u8>, config_json: &str) -> JsValue {
    let start = js_sys::Date::now() as u64;
    let reply = match run_analyze_cram(cram_bytes, crai_bytes, config_json) {
        Ok(mut ok) => {
            ok.elapsed_ms = (js_sys::Date::now() as u64).saturating_sub(start);
            AnalyzeReply::Ok(ok)
        }
        Err(message) => AnalyzeReply::Error { message },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

/// Run analysis on sparse byte ranges from File.slice() - BAM path.
/// `data`: concatenated bytes from all File.slice() calls.
/// `ranges_json`: `[{"file_offset": N, "length": N}, ...]` matching the data.
#[wasm_bindgen]
pub fn analyze_bam_sparse(
    data: Vec<u8>,
    ranges_json: &str,
    bai_bytes: Vec<u8>,
    config_json: &str,
) -> JsValue {
    let start = js_sys::Date::now() as u64;
    let reply = match run_analyze_sparse(data, ranges_json, bai_bytes, config_json, false) {
        Ok(mut ok) => {
            ok.elapsed_ms = (js_sys::Date::now() as u64).saturating_sub(start);
            AnalyzeReply::Ok(ok)
        }
        Err(message) => AnalyzeReply::Error { message },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

/// Run analysis on sparse byte ranges from File.slice() - CRAM path.
/// `ref_chr17` and `ref_chr2`: raw reference sequence bytes (no FASTA header).
#[wasm_bindgen]
pub fn analyze_cram_sparse(
    data: Vec<u8>,
    ranges_json: &str,
    crai_bytes: Vec<u8>,
    ref_chr17: Vec<u8>,
    ref_chr2: Vec<u8>,
    config_json: &str,
) -> JsValue {
    let start = js_sys::Date::now() as u64;
    let ref_seqs = vec![
        ("chr17".to_string(), ref_chr17),
        ("chr2".to_string(), ref_chr2),
    ];
    let reply = match run_analyze_sparse_cram(data, ranges_json, crai_bytes, ref_seqs, config_json)
    {
        Ok(mut ok) => {
            ok.elapsed_ms = (js_sys::Date::now() as u64).saturating_sub(start);
            AnalyzeReply::Ok(ok)
        }
        Err(message) => AnalyzeReply::Error { message },
    };
    serde_wasm_bindgen::to_value(&reply).unwrap_or(JsValue::NULL)
}

fn run_analyze_sparse(
    data: Vec<u8>,
    ranges_json: &str,
    index_bytes: Vec<u8>,
    config_json: &str,
    is_cram: bool,
) -> Result<AnalyzeOk, String> {
    let config = parse_config(config_json)?;

    // Parse the range descriptors - DON'T copy the data buffer.
    // Use OffsetCursor::from_concatenated for zero-copy construction.
    let range_descs: Vec<SparseRange> =
        serde_json::from_str(ranges_json).map_err(|e| format!("invalid ranges JSON: {e}"))?;
    let descs: Vec<(u64, u64)> = range_descs
        .iter()
        .map(|r| (r.file_offset, r.length))
        .collect();

    use nds_core::sources::offset_cursor::OffsetCursor;
    let cursor = OffsetCursor::from_concatenated(data, &descs);

    if is_cram {
        let mut source = CursorCramSource::from_cursor(cursor, index_bytes)
            .map_err(|e| format!("open CRAM failed: {e}"))?;
        run_analyze_source(&mut source, config)
    } else {
        let mut source = CursorBamSource::from_cursor(cursor, index_bytes)
            .map_err(|e| format!("open BAM failed: {e}"))?;
        run_analyze_source(&mut source, config)
    }
}

fn run_analyze_sparse_cram(
    data: Vec<u8>,
    ranges_json: &str,
    crai_bytes: Vec<u8>,
    ref_seqs: Vec<(String, Vec<u8>)>,
    config_json: &str,
) -> Result<AnalyzeOk, String> {
    let config = parse_config(config_json)?;
    let range_descs: Vec<SparseRange> =
        serde_json::from_str(ranges_json).map_err(|e| format!("invalid ranges JSON: {e}"))?;
    let descs: Vec<(u64, u64)> = range_descs
        .iter()
        .map(|r| (r.file_offset, r.length))
        .collect();

    use nds_core::sources::offset_cursor::OffsetCursor;
    let cursor = OffsetCursor::from_concatenated(data, &descs);
    let mut source = CursorCramSource::from_cursor_with_ref(cursor, crai_bytes, ref_seqs)
        .map_err(|e| format!("open CRAM failed: {e}"))?;
    run_analyze_source(&mut source, config)
}

#[derive(serde::Deserialize)]
struct SparseRange {
    file_offset: u64,
    length: u64,
}

fn run_analyze_cram(
    cram_bytes: Vec<u8>,
    crai_bytes: Vec<u8>,
    config_json: &str,
) -> Result<AnalyzeOk, String> {
    let config = parse_config(config_json)?;
    let mut source =
        CursorCramSource::new(cram_bytes, crai_bytes).map_err(|e| format!("open failed: {e}"))?;
    run_analyze_source(&mut source, config)
}

fn run_analyze(
    bam_bytes: Vec<u8>,
    bai_bytes: Vec<u8>,
    config_json: &str,
) -> Result<AnalyzeOk, String> {
    let config = parse_config(config_json)?;
    let mut source =
        CursorBamSource::new(bam_bytes, bai_bytes).map_err(|e| format!("open failed: {e}"))?;
    run_analyze_source(&mut source, config)
}

fn run_analyze_source<S: AlignmentSource>(
    source: &mut S,
    config: AnalysisConfig,
) -> Result<AnalyzeOk, String> {
    let header: AlignmentHeader = source
        .read_header()
        .map_err(|e| format!("read_header failed: {e}"))?;
    let chr17_md5 = header.chr17_md5().unwrap_or([0u8; 16]);
    let mut reference = InlineReferenceSource { chr17_md5 };

    let mut progress_events: Vec<ProgressDto> = Vec::new();
    let result = analyze_with_progress(source, &mut reference, &config, |event| {
        progress_events.push(event.into());
    })
    .map_err(|e| format!("analyze failed (exit {}): {}", e.exit_code(), e))?;

    let dup = result
        .boundary
        .windows
        .iter()
        .filter(|w| w.normalized_depth > config.duplication_threshold)
        .count();
    let del = result
        .boundary
        .windows
        .iter()
        .filter(|w| w.normalized_depth < config.deletion_threshold)
        .count();
    let runs: Vec<BoundaryRunDto> = result
        .boundary
        .runs
        .iter()
        .map(|run| BoundaryRunDto {
            start: run.start,
            end: run.end,
            length: run.length,
            direction: format!("{:?}", run.direction),
            mean_normalized_depth: run.mean_normalized_depth,
            bridged_gaps: run.bridged_gaps,
        })
        .collect();

    Ok(AnalyzeOk {
        depth: DepthDto {
            pmp22_mean_depth: result.depth.pmp22_mean_depth,
            autosomal_mean_depth: result.depth.autosomal_mean_depth,
            ratio: result.depth.ratio,
            estimated_cn: result.depth.estimated_cn,
        },
        boundary: BoundaryDto {
            window_count: result.boundary.windows.len(),
            duplicated_windows: dup,
            deleted_windows: del,
            runs,
        },
        interpretation: InterpretationDto {
            copy_number: result.interpretation.copy_number,
            confidence: format!("{:?}", result.interpretation.confidence),
            subtype_call: result.interpretation.subtype_call.map(|s| format!("{s:?}")),
            plain_language: result.interpretation.plain_language,
        },
        progress: progress_events,
        header_sequences: header.sequences.len(),
        elapsed_ms: 0,
    })
}

fn parse_config(json: &str) -> Result<AnalysisConfig, String> {
    if json.is_empty() || json == "{}" {
        return Ok(AnalysisConfig::default());
    }
    let overrides: ConfigOverrides =
        serde_json::from_str(json).map_err(|e| format!("invalid config JSON: {e}"))?;
    let def = AnalysisConfig::default();
    Ok(AnalysisConfig {
        chr2_control_region: Region::new(
            "chr2",
            overrides
                .chr2_control_start
                .unwrap_or(def.chr2_control_region.start),
            overrides
                .chr2_control_end
                .unwrap_or(def.chr2_control_region.end),
        ),
        boundary_scan_region: Region::new(
            "chr17",
            overrides
                .boundary_scan_start
                .unwrap_or(def.boundary_scan_region.start),
            overrides
                .boundary_scan_end
                .unwrap_or(def.boundary_scan_region.end),
        ),
        ..def
    })
}

// ---------------------------------------------------------------------------
// Shared types
// ---------------------------------------------------------------------------

struct InlineReferenceSource {
    chr17_md5: [u8; 16],
}

impl ReferenceSource for InlineReferenceSource {
    fn fetch_sequence(&mut self, _region: &Region) -> CoreResult<Vec<u8>> {
        Ok(Vec::new())
    }
    fn chr17_md5(&mut self) -> CoreResult<[u8; 16]> {
        Ok(self.chr17_md5)
    }
}

#[derive(serde::Deserialize)]
struct RegionInput {
    chrom: String,
    start: u32,
    end: u32,
}

#[derive(serde::Deserialize, Default)]
struct ConfigOverrides {
    chr2_control_start: Option<u32>,
    chr2_control_end: Option<u32>,
    boundary_scan_start: Option<u32>,
    boundary_scan_end: Option<u32>,
}

#[derive(Serialize)]
struct RangeDto {
    offset: u64,
    length: u64,
}

// --- Analysis DTOs ---

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum ProgressDto {
    HeaderVerified,
    WalkStart {
        phase: String,
        region: String,
    },
    WalkTick {
        phase: String,
        records_seen: u64,
        last_position: u32,
    },
    WalkDone {
        phase: String,
        records_seen: u64,
    },
    DepthRatioDone {
        pmp22_mean_depth: f64,
        autosomal_mean_depth: f64,
        ratio: f64,
        estimated_cn: u8,
    },
    BoundaryScanDone {
        windows: usize,
        longest_run_length: Option<u32>,
    },
    InterpretDone {
        confidence: String,
        subtype_call: Option<String>,
    },
    Unknown,
}

impl From<Progress> for ProgressDto {
    fn from(event: Progress) -> Self {
        match event {
            Progress::HeaderVerified => Self::HeaderVerified,
            Progress::WalkStart { phase, region } => Self::WalkStart {
                phase: phase.to_string(),
                region,
            },
            Progress::WalkTick {
                phase,
                records_seen,
                last_position,
            } => Self::WalkTick {
                phase: phase.to_string(),
                records_seen,
                last_position,
            },
            Progress::WalkDone {
                phase,
                records_seen,
            } => Self::WalkDone {
                phase: phase.to_string(),
                records_seen,
            },
            Progress::DepthRatioDone {
                pmp22_mean_depth,
                autosomal_mean_depth,
                ratio,
                estimated_cn,
            } => Self::DepthRatioDone {
                pmp22_mean_depth,
                autosomal_mean_depth,
                ratio,
                estimated_cn,
            },
            Progress::BoundaryScanDone {
                windows,
                longest_run_length,
            } => Self::BoundaryScanDone {
                windows,
                longest_run_length,
            },
            Progress::InterpretDone {
                confidence,
                subtype_call,
            } => Self::InterpretDone {
                confidence: format!("{confidence:?}"),
                subtype_call: subtype_call.map(|s| format!("{s:?}")),
            },
            _ => Self::Unknown,
        }
    }
}

#[derive(Serialize)]
struct DepthDto {
    pmp22_mean_depth: f64,
    autosomal_mean_depth: f64,
    ratio: f64,
    estimated_cn: u8,
}

#[derive(Serialize)]
struct BoundaryRunDto {
    start: u32,
    end: u32,
    length: u32,
    direction: String,
    mean_normalized_depth: f64,
    bridged_gaps: u32,
}

#[derive(Serialize)]
struct BoundaryDto {
    window_count: usize,
    duplicated_windows: usize,
    deleted_windows: usize,
    runs: Vec<BoundaryRunDto>,
}

#[derive(Serialize)]
struct InterpretationDto {
    copy_number: Option<u8>,
    confidence: String,
    subtype_call: Option<String>,
    plain_language: String,
}

#[derive(Serialize)]
struct AnalyzeOk {
    depth: DepthDto,
    boundary: BoundaryDto,
    interpretation: InterpretationDto,
    progress: Vec<ProgressDto>,
    header_sequences: usize,
    elapsed_ms: u64,
}

#[derive(Serialize)]
#[serde(tag = "status")]
enum AnalyzeReply {
    #[serde(rename = "ok")]
    Ok(AnalyzeOk),
    #[serde(rename = "error")]
    Error { message: String },
}
