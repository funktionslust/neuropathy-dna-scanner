//! Generate test BAM fixtures from the ClinVar catalog.
//!
//! Reads catalog.json, picks representative variants, generates BAM files
//! with those variants spiked in at the correct positions. The output is
//! a set of small BAM files + a JSON manifest listing what was spiked.
//!
//! Usage:
//!   build_test_fixtures --catalog catalog.json --output-dir fixtures/

use clap::Parser;
use nds_core::catalog::VariantCatalog;
use nds_tools::synth::{self, CoverageSpec, VariantSpike};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "build_test_fixtures",
    about = "Generate test CRAMs from ClinVar catalog"
)]
struct Cli {
    /// Path to catalog.json (from build_catalog).
    #[arg(long)]
    catalog: PathBuf,

    /// Reference FASTA for CRAM generation.
    #[arg(long)]
    reference: PathBuf,

    /// Output directory for fixtures.
    #[arg(long, default_value = "fixtures")]
    output_dir: PathBuf,
}

/// A test scenario: what variants to spike and what SV profile to create.
#[derive(Debug, serde::Serialize)]
struct FixtureManifest {
    fixtures: Vec<FixtureEntry>,
}

#[derive(Debug, serde::Serialize)]
struct FixtureEntry {
    name: String,
    description: String,
    cram_file: String,
    /// SV profile: "cn2" (normal), "cn3_cmt1a", "cn1_hnpp", etc.
    sv_profile: String,
    /// Spiked SNV variants with expected calls.
    spiked_variants: Vec<SpikedVariant>,
}

#[derive(Debug, Clone, serde::Serialize)]
struct SpikedVariant {
    gene: String,
    chrom: String,
    pos: u32,
    ref_allele: String,
    alt_allele: String,
    clinvar_id: u32,
    expected_call: String, // "heterozygous", "homozygous", "reference_only"
}

fn main() {
    let cli = Cli::parse();

    let catalog_data = std::fs::read_to_string(&cli.catalog).expect("failed to read catalog");
    let catalog: VariantCatalog =
        serde_json::from_str(&catalog_data).expect("failed to parse catalog");

    eprintln!(
        "Loaded catalog: {} variants across {} genes",
        catalog.meta.variant_count, catalog.meta.gene_count
    );

    std::fs::create_dir_all(&cli.output_dir).expect("failed to create output dir");

    // Pick representative variants from the catalog for each CMT subtype.
    // Strategy: for each key gene, pick the variant with the most review stars.
    let cmt_genes = [
        ("PMP22", "CMT1E / HNPP point mutation"),
        ("GJB1", "CMTX1"),
        ("MPZ", "CMT1B"),
        ("MFN2", "CMT2A"),
        ("GDAP1", "CMT4A"),
        ("SH3TC2", "CMT4C"),
        ("NEFL", "CMT2E"),
        ("EGR2", "CMT1D"),
    ];

    let mut het_variants: Vec<(String, SpikedVariant, VariantSpike)> = Vec::new();

    for &(gene, subtype) in &cmt_genes {
        // Find the best SNV variant for this gene (highest stars, SNV only)
        let best = catalog
            .variants
            .iter()
            .filter(|v| v.gene == gene)
            .filter(|v| v.ref_allele.len() == 1 && v.alt_allele.len() == 1)
            .max_by_key(|v| v.review_stars);

        if let Some(v) = best {
            eprintln!(
                "  {gene} ({subtype}): {chr}:{pos} {ref_a}>{alt_a} ClinVar:{id} ({stars} stars)",
                chr = v.chrom,
                pos = v.pos,
                ref_a = v.ref_allele,
                alt_a = v.alt_allele,
                id = v.clinvar_id,
                stars = v.review_stars
            );

            let spike = VariantSpike::het(
                format!("chr{}", v.chrom),
                v.pos,
                v.ref_allele.as_bytes()[0],
                v.alt_allele.as_bytes()[0],
            );
            let spiked = SpikedVariant {
                gene: gene.to_string(),
                chrom: v.chrom.clone(),
                pos: v.pos,
                ref_allele: v.ref_allele.clone(),
                alt_allele: v.alt_allele.clone(),
                clinvar_id: v.clinvar_id,
                expected_call: "heterozygous".to_string(),
            };
            het_variants.push((subtype.to_string(), spiked, spike));
        } else {
            eprintln!("  {gene} ({subtype}): no SNV variants in catalog");
        }
    }

    let mut fixtures: Vec<FixtureEntry> = Vec::new();

    // Fixture 1: CN=2 normal + all CMT point mutations (heterozygous)
    // Coverage specs are derived from the actual selected variant positions
    // so every spike is guaranteed to land inside a covered region.
    {
        let name = "cn2_with_snv_het";
        let mut specs = normal_coverage_specs();
        // Add coverage around each selected variant position (±2 kb)
        for (_, sv, _) in &het_variants {
            let chrom = format!("chr{}", sv.chrom);
            // chr17 is already covered by the base spec.
            if chrom == "chr17" {
                continue;
            }
            specs.push(CoverageSpec {
                chrom,
                start: sv.pos.saturating_sub(2000),
                end: sv.pos + 2000,
                depth: 30,
            });
        }
        let spikes: Vec<VariantSpike> = het_variants.iter().map(|(_, _, s)| s.clone()).collect();
        let spiked: Vec<SpikedVariant> = het_variants.iter().map(|(_, sv, _)| sv.clone()).collect();

        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &spikes, 42);
        eprintln!(
            "Generated {name}.cram with {} spiked variants",
            spikes.len()
        );

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Normal CN=2 with heterozygous CMT point mutations spiked in".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: spiked,
        });
    }

    // Fixture 2: CN=3 CMT1A duplication (no SNV spikes)
    // Tests SV detection without SNV findings.
    {
        let name = "cn3_cmt1a_sv_only";
        let specs = cmt1a_coverage_specs();
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 43);
        eprintln!("Generated {name}.cram (CMT1A duplication, no SNV)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "CMT1A duplication (CN=3) without point mutations".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn3_cmt1a".to_string(),
            spiked_variants: vec![],
        });
    }

    // Fixture 3: CN=2 normal, no variants (negative control)
    {
        let name = "cn2_negative_control";
        let specs = normal_coverage_specs();
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 44);
        eprintln!("Generated {name}.cram (negative control)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Normal CN=2, no variants spiked (negative control)".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![],
        });
    }

    // Fixture 4: CN=1 HNPP deletion (no SNV spikes)
    {
        let name = "cn1_hnpp_sv_only";
        let specs = hnpp_coverage_specs();
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 45);
        eprintln!("Generated {name}.cram (HNPP deletion)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "HNPP deletion (CN=1) without point mutations".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn1_hnpp".to_string(),
            spiked_variants: vec![],
        });
    }

    // Regression fixtures covering the most frequent false-positive paths.

    // Strand-bias + low base-quality must NOT produce a heterozygous call.
    // Placed in the chr17 PMP22 region to avoid cross-chrom CRAM issues.
    {
        let name = "cn2_pnpla6_strand_bias";
        let specs = normal_coverage_specs();
        // Spike T>G with strand bias + low BQ at a chr17 position
        // (using PMP22 region to avoid cross-chrom CRAM issues)
        let spikes = vec![VariantSpike {
            chrom: "chr17".into(),
            pos: 15_240_000,
            ref_base: b'T',
            alt_base: b'G',
            allele_fraction: 0.3,
            strand_override: Some(synth::Strand::Reverse),
            alt_base_quality: Some(10), // below min_base_qual=20
            alt_mapq: None,
        }];
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &spikes, 50);
        eprintln!("Generated {name} (PNPLA6 strand bias regression)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description:
                "PNPLA6 chr19:7555369 T>G with strand bias + low BQ. Must NOT call heterozygous."
                    .to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![SpikedVariant {
                gene: "strand_bias_test".to_string(),
                chrom: "17".to_string(),
                pos: 15_240_000,
                ref_allele: "T".to_string(),
                alt_allele: "G".to_string(),
                clinvar_id: 2683056,
                expected_call: "no_call".to_string(), // must NOT be heterozygous
            }],
        });
    }

    // Segdup artifact - elevated depth at SMS-REP middle
    {
        let name = "cn2_segdup_artifact";
        let mut specs = normal_coverage_specs();
        // Add artificial 2x elevation at the known segdup
        specs.push(CoverageSpec {
            chrom: "chr17".into(),
            start: 18_597_985,
            end: 18_622_233,
            depth: 30, // adds to existing 30x = 60x total = 2x elevation
        });
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 51);
        eprintln!("Generated {name} (segdup artifact regression)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "SMS-REP middle 2x depth artifact. Must NOT call as duplication."
                .to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![],
        });
    }

    // Subthreshold run - tiny elevated region
    {
        let name = "cn2_subthreshold_run";
        let mut specs = normal_coverage_specs();
        // 12 kb elevated region inside the boundary scan area
        specs.push(CoverageSpec {
            chrom: "chr17".into(),
            start: 15_100_000,
            end: 15_112_000,
            depth: 15, // adds to 30x = 45x total = 1.5x
        });
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 52);
        eprintln!("Generated {name} (subthreshold run regression)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "12 kb elevated region. Must NOT call as duplication (below min_reportable_run_length).".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![],
        });
    }

    // Multi-alt at same position
    // Find a ClinVar position with two different alt alleles
    {
        let name = "cn2_multi_alt_same_pos";
        // Find two variants at the same position in the catalog
        let mut multi_pos: Option<(u32, &str, u8, u8, u8)> = None; // (pos, chrom, ref, alt1, alt2)
        let mut seen: std::collections::HashMap<(String, u32), Vec<u8>> =
            std::collections::HashMap::new();
        for v in &catalog.variants {
            // Only chr17 to avoid cross-chrom CRAM issues in fixtures
            if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 && v.chrom == "17" {
                let key = (v.chrom.clone(), v.pos);
                let alt = v.alt_allele.as_bytes()[0];
                seen.entry(key).or_default().push(alt);
            }
        }
        for ((chrom, pos), alts) in &seen {
            if alts.len() >= 2 {
                let v = catalog
                    .variants
                    .iter()
                    .find(|v| v.chrom == *chrom && v.pos == *pos)
                    .unwrap();
                multi_pos = Some((*pos, &v.chrom, v.ref_allele.as_bytes()[0], alts[0], alts[1]));
                break;
            }
        }

        if let Some((pos, chrom, ref_base, alt1, _alt2)) = multi_pos {
            let mut specs = normal_coverage_specs();
            specs.push(CoverageSpec {
                chrom: format!("chr{chrom}"),
                start: pos.saturating_sub(500),
                end: pos + 500,
                depth: 30,
            });
            // Spike only alt1, not alt2
            let spikes = vec![VariantSpike::het(
                format!("chr{chrom}"),
                pos,
                ref_base,
                alt1,
            )];
            let path = cli.output_dir.join(format!("{name}.cram"));
            generate_fixture(&path, &cli.reference, &specs, &spikes, 53);
            eprintln!("Generated {name} (multi-alt at {chrom}:{pos})");

            fixtures.push(FixtureEntry {
                name: name.to_string(),
                description: format!("Two ClinVar variants at {chrom}:{pos}. Only alt1 spiked. Alt2 must NOT be called."),
                cram_file: format!("{name}.cram"),
                sv_profile: "cn2".to_string(),
                spiked_variants: vec![SpikedVariant {
                    gene: "multi_alt_test".to_string(),
                    chrom: chrom.to_string(),
                    pos,
                    ref_allele: String::from(ref_base as char),
                    alt_allele: String::from(alt1 as char),
                    clinvar_id: 0,
                    expected_call: "heterozygous".to_string(),
                }],
            });
        } else {
            eprintln!("  skipped {name}: no multi-alt position found in catalog");
        }
    }

    // Complex dup-normal-dup gap
    {
        let name = "cn2_complex_dup_gap";
        let specs = vec![
            CoverageSpec {
                chrom: "chr2".into(),
                start: 50_000_000,
                end: 50_500_000,
                depth: 30,
            },
            CoverageSpec {
                chrom: "chr17".into(),
                start: 13_500_000,
                end: 14_200_000,
                depth: 30,
            },
            // First dup segment
            CoverageSpec {
                chrom: "chr17".into(),
                start: 14_200_000,
                end: 14_700_000,
                depth: 45,
            },
            // Normal gap
            CoverageSpec {
                chrom: "chr17".into(),
                start: 14_700_000,
                end: 14_900_000,
                depth: 30,
            },
            // Second dup segment
            CoverageSpec {
                chrom: "chr17".into(),
                start: 14_900_000,
                end: 15_400_000,
                depth: 45,
            },
            CoverageSpec {
                chrom: "chr17".into(),
                start: 15_400_000,
                end: 21_000_000,
                depth: 30,
            },
        ];
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 54);
        eprintln!("Generated {name} (complex dup with gap)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Two 500kb dup segments separated by 200kb normal. Should NOT classify as standard CMT1A.".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn3_atypical".to_string(),
            spiked_variants: vec![],
        });
    }

    // Sex-aware chrX fixtures.

    // Male chrX hemizygous reference (no variant, adequate coverage)
    {
        let name = "cn2_male_chrx_hemi_ref";
        let mut specs = normal_coverage_specs();
        // chrX GJB1 region at half depth (male = 1 copy)
        specs.push(CoverageSpec {
            chrom: "chrX".into(),
            start: 71_221_000,
            end: 71_226_000,
            depth: 15, // ~half of autosomal 30x
        });
        // chrY at half depth (confirms male)
        specs.push(CoverageSpec {
            chrom: "chrY".into(),
            start: 2_786_000,
            end: 2_790_000,
            depth: 15,
        });
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &[], 55);
        eprintln!("Generated {name} (male chrX reference)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Male sample, chrX GJB1 region at hemizygous depth. No variant spiked."
                .to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![],
        });
    }

    // Male chrX hemizygous alt (GJB1 variant present)
    {
        let name = "cn2_male_chrx_hemi_alt";
        let mut specs = normal_coverage_specs();
        specs.push(CoverageSpec {
            chrom: "chrX".into(),
            start: 71_221_000,
            end: 71_226_000,
            depth: 15,
        });
        specs.push(CoverageSpec {
            chrom: "chrY".into(),
            start: 2_786_000,
            end: 2_790_000,
            depth: 15,
        });
        // Spike GJB1 C>T at chrX:71223249 with high AF (hemizygous = ~100%)
        let spikes = vec![VariantSpike {
            chrom: "chrX".into(),
            pos: 71_223_249,
            ref_base: b'C',
            alt_base: b'T',
            allele_fraction: 0.95, // hemizygous: nearly all reads show alt
            strand_override: None,
            alt_base_quality: None,
            alt_mapq: None,
        }];
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &spikes, 56);
        eprintln!("Generated {name} (male chrX hemizygous alt)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Male sample with hemizygous GJB1 C>T variant on chrX.".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![SpikedVariant {
                gene: "GJB1".to_string(),
                chrom: "X".to_string(),
                pos: 71_223_249,
                ref_allele: "C".to_string(),
                alt_allele: "T".to_string(),
                clinvar_id: 217166,
                expected_call: "hemizygous".to_string(),
            }],
        });
    }

    // Female chrX heterozygous (GJB1 variant, diploid depth)
    {
        let name = "cn2_female_chrx_het";
        let mut specs = normal_coverage_specs();
        // chrX at diploid depth (female = 2 copies)
        specs.push(CoverageSpec {
            chrom: "chrX".into(),
            start: 71_221_000,
            end: 71_226_000,
            depth: 30,
        });
        // NO chrY (female)
        // Spike GJB1 C>T at chrX:71223249 as het (~50% AF)
        let spikes = vec![VariantSpike::het(
            String::from("chrX"),
            71_223_249,
            b'C',
            b'T',
        )];
        let path = cli.output_dir.join(format!("{name}.cram"));
        generate_fixture(&path, &cli.reference, &specs, &spikes, 57);
        eprintln!("Generated {name} (female chrX het)");

        fixtures.push(FixtureEntry {
            name: name.to_string(),
            description: "Female sample with heterozygous GJB1 C>T variant on chrX.".to_string(),
            cram_file: format!("{name}.cram"),
            sv_profile: "cn2".to_string(),
            spiked_variants: vec![SpikedVariant {
                gene: "GJB1".to_string(),
                chrom: "X".to_string(),
                pos: 71_223_249,
                ref_allele: "C".to_string(),
                alt_allele: "T".to_string(),
                clinvar_id: 217166,
                expected_call: "heterozygous".to_string(),
            }],
        });
    }

    // Write manifest
    let manifest = FixtureManifest { fixtures };
    let manifest_path = cli.output_dir.join("manifest.json");
    let manifest_json = serde_json::to_string_pretty(&manifest).expect("serialize manifest");
    std::fs::write(&manifest_path, &manifest_json).expect("write manifest");
    eprintln!("\nManifest written to {}", manifest_path.display());
}

fn generate_fixture(
    cram_path: &std::path::Path,
    reference: &std::path::Path,
    specs: &[CoverageSpec],
    spikes: &[VariantSpike],
    seed: u64,
) {
    let header = synth::build_grch38_header();

    // Generate reads with actual reference bases (CRAM-compatible)
    let reads = synth::generate_reads_with_ref(&header, specs, spikes, reference, seed);

    // Write CRAM directly using noodles
    synth::write_cram(cram_path, &header, &reads, reference).expect("write CRAM failed");

    // Index CRAM -> CRAI using samtools
    let status = std::process::Command::new("samtools")
        .args(["index", cram_path.to_str().unwrap()])
        .status()
        .expect("samtools not found");
    if !status.success() {
        eprintln!("  ERROR: samtools index failed for {}", cram_path.display());
    }
}

/// Base SV coverage: chr17 + chr2 only. Safe for all fixtures because the
/// browser only loads reference for these two chromosomes during SV analysis.
/// The broad chr17 region covers PMP22 + boundary scan + RAI1, so fixtures
/// must not add a separate PMP22 spec (would double-count reads).
fn normal_coverage_specs() -> Vec<CoverageSpec> {
    vec![
        CoverageSpec {
            chrom: "chr17".into(),
            start: 13_500_000,
            end: 21_000_000,
            depth: 30,
        },
        CoverageSpec {
            chrom: "chr2".into(),
            start: 50_000_000,
            end: 50_500_000,
            depth: 30,
        },
    ]
}

fn cmt1a_coverage_specs() -> Vec<CoverageSpec> {
    // CMT1A: 1.4 Mb duplication between CMT1A-REPs (GRCh38: ~14.17M - ~15.59M)
    // Depth 45 = 1.5x baseline (CN=3). PMP22 is inside the dup region.
    let dup_start = 14_170_000;
    let dup_end = 15_590_000;
    vec![
        CoverageSpec {
            chrom: "chr2".into(),
            start: 50_000_000,
            end: 50_500_000,
            depth: 30,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: 13_500_000,
            end: dup_start,
            depth: 30,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: dup_start,
            end: dup_end,
            depth: 45,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: dup_end,
            end: 21_000_000,
            depth: 30,
        },
    ]
}

fn hnpp_coverage_specs() -> Vec<CoverageSpec> {
    // HNPP: reciprocal deletion of the same 1.4 Mb CMT1A region
    // Depth 15 = 0.5x baseline (CN=1). PMP22 is inside the del region
    // so it inherits 15x automatically - no separate PMP22 spec needed.
    let del_start = 14_170_000;
    let del_end = 15_590_000;
    vec![
        CoverageSpec {
            chrom: "chr2".into(),
            start: 50_000_000,
            end: 50_500_000,
            depth: 30,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: 13_500_000,
            end: del_start,
            depth: 30,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: del_start,
            end: del_end,
            depth: 15,
        },
        CoverageSpec {
            chrom: "chr17".into(),
            start: del_end,
            end: 21_000_000,
            depth: 30,
        },
    ]
}
