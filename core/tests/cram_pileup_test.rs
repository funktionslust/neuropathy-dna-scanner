//! End-to-end CRAM → pileup → classification tests.
//!
//! These tests load CRAM fixtures, run them through the full pileup pipeline
//! (core::pileup_cram::pileup_cram), and assert per-variant calls.
//! This guards against regressions in:
//!   - base_at_pos_cram (CIGAR walking, base extraction)
//!   - BQ filter (min_base_qual application)
//!   - Strand tracking (forward/reverse counting)
//!   - Alt-allele matching (per-variant keying)
//!   - Chromosome name resolution (ClinVar "17" → CRAM "chr17")
//!
//! Requires GRCh38 reference — skipped if not present.

use std::fs::File;
use std::io::{BufReader, Cursor, Read as _};
use std::path::PathBuf;

use nds_core::catalog::{CatalogMeta, ClinSig, Inheritance, PathogenicVariant, VariantCatalog};
use nds_core::pileup::{PileupConfig, VariantCall};
use nds_core::pileup_cram;
use nds_core::Sex;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fixture_path(name: &str) -> PathBuf {
    repo_root().join("fixtures").join(name)
}

fn reference_path() -> PathBuf {
    repo_root().join("ref/GRCh38_no_alt_plus_hs38d1.fna")
}

fn reference_available() -> bool {
    reference_path().exists()
}

fn load_bytes(path: &PathBuf) -> Vec<u8> {
    let mut bytes = Vec::new();
    File::open(path)
        .unwrap_or_else(|e| panic!("open {}: {e}", path.display()))
        .read_to_end(&mut bytes)
        .unwrap();
    bytes
}

/// Build a minimal catalog with specific test variants.
fn test_catalog(variants: Vec<PathogenicVariant>) -> VariantCatalog {
    let count = variants.len();
    VariantCatalog {
        meta: CatalogMeta {
            generated_date: "2024-01-01".into(),
            clinvar_date: "test".into(),
            panelapp_panel: "test".into(),
            gene_count: 1,
            variant_count: count,
            min_review_stars: 0,
            excluded_genes: vec![],
        },
        variants,
    }
}

fn test_variant(chrom: &str, pos: u32, ref_a: &str, alt_a: &str) -> PathogenicVariant {
    PathogenicVariant {
        gene: "TEST".into(),
        chrom: chrom.into(),
        pos,
        ref_allele: ref_a.into(),
        alt_allele: alt_a.into(),
        clinvar_id: 0,
        significance: ClinSig::Pathogenic,
        review_stars: 2,
        condition: "test".into(),
        inheritance: Inheritance::Unknown,
        consequence: "missense".into(),
    }
}

/// Open a CRAM fixture with reference and run pileup against a catalog.
fn run_pileup_on_fixture(
    cram_name: &str,
    catalog: &VariantCatalog,
) -> Option<nds_core::pileup::PileupResult> {
    run_pileup_on_fixture_with_sex(cram_name, catalog, Sex::Unknown)
}

fn run_pileup_on_fixture_with_sex(
    cram_name: &str,
    catalog: &VariantCatalog,
    sex: Sex,
) -> Option<nds_core::pileup::PileupResult> {
    if !reference_available() {
        return None;
    }

    let cram_path = fixture_path(cram_name);
    let crai_path = fixture_path(&format!("{cram_name}.crai"));
    if !cram_path.exists() {
        panic!("Fixture not found: {}", cram_path.display());
    }

    let cram_bytes = load_bytes(&cram_path);
    let crai_bytes = load_bytes(&crai_path);

    // Build reference sequences for chromosomes in the catalog
    let ref_path = reference_path();
    let mut fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
        .build_from_path(&ref_path)
        .unwrap();

    // Collect unique chromosomes from catalog
    let chroms: std::collections::HashSet<String> = catalog
        .variants
        .iter()
        .map(|v| {
            if v.chrom.starts_with("chr") {
                v.chrom.clone()
            } else {
                format!("chr{}", v.chrom)
            }
        })
        .collect();

    // GRCh38 chromosome lengths (for the ones we need)
    let chrom_lengths: std::collections::HashMap<&str, usize> = [
        ("chr1", 248956422),
        ("chr2", 242193529),
        ("chr5", 181538259),
        ("chr8", 145138636),
        ("chr10", 133797422),
        ("chr17", 83257441),
        ("chrX", 156040895),
        ("chrY", 57227415),
    ]
    .into_iter()
    .collect();

    let mut ref_seqs: Vec<(String, Vec<u8>)> = Vec::new();
    for chrom in &chroms {
        let len = match chrom_lengths.get(chrom.as_str()) {
            Some(&l) => l,
            None => continue,
        };
        let region = noodles_core::Region::new(
            chrom.as_str(),
            noodles_core::Position::try_from(1usize).unwrap()
                ..=noodles_core::Position::try_from(len).unwrap(),
        );
        let record = fasta_reader.query(&region).unwrap();
        ref_seqs.push((chrom.clone(), record.sequence().as_ref().to_vec()));
    }

    // Build CRAM reader with reference
    let crai_index = noodles_cram::crai::io::Reader::new(Cursor::new(crai_bytes))
        .read_index()
        .unwrap();

    let fasta_records: Vec<noodles_fasta::Record> = ref_seqs
        .into_iter()
        .map(|(name, seq_bytes)| {
            let def = noodles_fasta::record::Definition::new(name, None);
            let seq = noodles_fasta::record::Sequence::from(seq_bytes);
            noodles_fasta::Record::new(def, seq)
        })
        .collect();
    let repo = noodles_fasta::repository::Repository::new(fasta_records);

    let mut reader = noodles_cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_reader(BufReader::new(Cursor::new(cram_bytes)));
    let header = reader.read_header().unwrap();

    let config = PileupConfig::default();
    let result =
        pileup_cram::pileup_cram(&mut reader, &header, &crai_index, catalog, &config, sex).unwrap();
    Some(result)
}

// =========================================================================
// cn2_with_snv_het: all 8 spiked variants must be called heterozygous
// =========================================================================

#[test]
fn pileup_cn2_with_snv_het_finds_all_spikes() {
    // Read the manifest to get the actual spiked variants
    let manifest_path = fixture_path("manifest.json");
    if !manifest_path.exists() {
        eprintln!("SKIP: fixtures not generated");
        return;
    }
    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest_path).unwrap()).unwrap();
    let fixture = manifest["fixtures"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "cn2_with_snv_het")
        .expect("cn2_with_snv_het not in manifest");

    let spikes = fixture["spiked_variants"].as_array().unwrap();
    let variants: Vec<PathogenicVariant> = spikes
        .iter()
        .map(|s| {
            test_variant(
                s["chrom"].as_str().unwrap(),
                s["pos"].as_u64().unwrap() as u32,
                s["ref_allele"].as_str().unwrap(),
                s["alt_allele"].as_str().unwrap(),
            )
        })
        .collect();

    let catalog = test_catalog(variants);
    let Some(result) = run_pileup_on_fixture("cn2_with_snv_het.cram", &catalog) else {
        eprintln!("SKIP: reference not available");
        return;
    };

    // All 8 spikes should be positive findings (heterozygous)
    assert_eq!(
        result.positive_findings.len(),
        spikes.len(),
        "expected {} positive findings, got {}. Findings: {:?}",
        spikes.len(),
        result.positive_findings.len(),
        result
            .positive_findings
            .iter()
            .map(|f| format!("{}:{} {:?}", f.chrom, f.pos, f.call))
            .collect::<Vec<_>>()
    );

    for finding in &result.positive_findings {
        assert!(
            matches!(finding.call, VariantCall::Heterozygous { .. }),
            "expected Heterozygous at {}:{}, got {:?}",
            finding.chrom,
            finding.pos,
            finding.call
        );
    }
}

// =========================================================================
// cn2_pnpla6_strand_bias: must NOT call het (BQ filter + strand bias)
// =========================================================================

#[test]
fn pileup_pnpla6_strand_bias_not_called() {
    // The spike is at chr17:15240000 T>G with BQ=10 + reverse-only
    let catalog = test_catalog(vec![test_variant("17", 15_240_000, "T", "G")]);

    let Some(result) = run_pileup_on_fixture("cn2_pnpla6_strand_bias.cram", &catalog) else {
        eprintln!("SKIP: reference not available");
        return;
    };

    // Must NOT be a positive finding
    assert_eq!(
        result.positive_findings.len(),
        0,
        "PNPLA6 strand bias spike must NOT be called positive. Got: {:?}",
        result
            .positive_findings
            .iter()
            .map(|f| format!("{}:{} {:?}", f.chrom, f.pos, f.call))
            .collect::<Vec<_>>()
    );
}

// =========================================================================
// cn2_multi_alt_same_pos: only the spiked alt should be called
// =========================================================================

#[test]
fn pileup_multi_alt_only_spiked_alt_called() {
    // Read manifest to get the actual position and alleles
    let manifest_path = fixture_path("manifest.json");
    if !manifest_path.exists() {
        eprintln!("SKIP: fixtures not generated");
        return;
    }
    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&manifest_path).unwrap()).unwrap();
    let fixture = manifest["fixtures"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "cn2_multi_alt_same_pos")
        .expect("cn2_multi_alt_same_pos not in manifest");

    let spike = &fixture["spiked_variants"][0];
    let pos = spike["pos"].as_u64().unwrap() as u32;
    let ref_a = spike["ref_allele"].as_str().unwrap();
    let spiked_alt = spike["alt_allele"].as_str().unwrap();

    // Create catalog with the spiked alt AND a different alt at same position.
    // Pick an alt that won't appear in real reads. Since ref is known,
    // use a base that differs from both ref and spiked alt.
    let bases = ["A", "T", "C", "G"];
    let other_alt = bases
        .iter()
        .find(|&&b| b != ref_a && b != spiked_alt)
        .unwrap();

    let catalog = test_catalog(vec![
        test_variant("17", pos, ref_a, spiked_alt),
        test_variant("17", pos, ref_a, other_alt),
    ]);

    let Some(result) = run_pileup_on_fixture("cn2_multi_alt_same_pos.cram", &catalog) else {
        eprintln!("SKIP: reference not available");
        return;
    };

    // The spiked alt must be positive
    let spiked_finding = result
        .positive_findings
        .iter()
        .find(|f| f.alt_allele == spiked_alt);
    assert!(
        spiked_finding.is_some(),
        "spiked alt {spiked_alt} at {pos} must be called positive"
    );

    // The other alt must NOT be positive
    let other_finding = result
        .positive_findings
        .iter()
        .find(|f| f.alt_allele == *other_alt);
    assert!(
        other_finding.is_none(),
        "unspiked alt {other_alt} at {pos} must NOT be called positive, but got {:?}",
        other_finding.map(|f| &f.call)
    );
}

// =========================================================================
// cn2_negative_control: no variants should be called at any position
// =========================================================================

#[test]
fn pileup_negative_control_finds_nothing() {
    // Check a few known ClinVar positions — none should have alt alleles
    let catalog = test_catalog(vec![
        test_variant("17", 15_260_693, "T", "C"), // PMP22
        test_variant("17", 15_240_000, "T", "G"), // arbitrary
    ]);

    let Some(result) = run_pileup_on_fixture("cn2_negative_control.cram", &catalog) else {
        eprintln!("SKIP: reference not available");
        return;
    };

    assert_eq!(
        result.positive_findings.len(),
        0,
        "negative control should have no positive findings"
    );
}

// =========================================================================
// Sex-aware chrX fixture tests (end-to-end through CRAM pileup path)
// =========================================================================

#[test]
fn pileup_male_chrx_hemizygous_alt_detected() {
    // Male fixture with GJB1 C>T at chrX:71223249, 95% AF
    let catalog = test_catalog(vec![test_variant("X", 71_223_249, "C", "T")]);

    let Some(result) =
        run_pileup_on_fixture_with_sex("cn2_male_chrx_hemi_alt.cram", &catalog, Sex::Male)
    else {
        eprintln!("SKIP: reference not available");
        return;
    };

    assert_eq!(
        result.positive_findings.len(),
        1,
        "male hemizygous alt should produce 1 positive finding"
    );
    assert!(
        matches!(
            result.positive_findings[0].call,
            VariantCall::Hemizygous { .. }
        ),
        "expected Hemizygous call for male chrX, got {:?}",
        result.positive_findings[0].call
    );
}

#[test]
fn pileup_male_chrx_reference_not_called() {
    // Male fixture with no variant spiked at GJB1 position
    let catalog = test_catalog(vec![test_variant("X", 71_223_249, "C", "T")]);

    let Some(result) =
        run_pileup_on_fixture_with_sex("cn2_male_chrx_hemi_ref.cram", &catalog, Sex::Male)
    else {
        eprintln!("SKIP: reference not available");
        return;
    };

    assert_eq!(
        result.positive_findings.len(),
        0,
        "male chrX reference-only should have no positive findings"
    );
}

#[test]
fn pileup_female_chrx_het_detected() {
    // Female fixture with GJB1 C>T at chrX:71223249, ~50% AF (diploid het)
    let catalog = test_catalog(vec![test_variant("X", 71_223_249, "C", "T")]);

    let Some(result) =
        run_pileup_on_fixture_with_sex("cn2_female_chrx_het.cram", &catalog, Sex::Female)
    else {
        eprintln!("SKIP: reference not available");
        return;
    };

    assert_eq!(
        result.positive_findings.len(),
        1,
        "female diploid het should produce 1 positive finding"
    );
    assert!(
        matches!(
            result.positive_findings[0].call,
            VariantCall::Heterozygous { .. }
        ),
        "expected Heterozygous call for female chrX, got {:?}",
        result.positive_findings[0].call
    );
}
