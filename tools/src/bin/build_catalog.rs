//! Builds the pathogenic variant catalog from ClinVar + PanelApp data.
//!
//! Usage:
//!   build_catalog --clinvar clinvar.vcf.gz --panelapp panel.json --output catalog.json
//!
//! Or fetch PanelApp live:
//!   build_catalog --clinvar clinvar.vcf.gz --panelapp-id 85 --output catalog.json
//!
//! The output is a JSON file that can be embedded in the core crate or
//! loaded at runtime. Every variant is traceable to its ClinVar ID.

use clap::Parser;
use flate2::read::MultiGzDecoder;
use nds_core::catalog::{
    CatalogMeta, ClinSig, ExcludedGene, ExclusionReason, Inheritance, PathogenicVariant,
    VariantCatalog,
};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "build_catalog",
    about = "Build CMT variant catalog from ClinVar + PanelApp"
)]
struct Cli {
    /// Path to clinvar.vcf.gz (download from NCBI FTP).
    #[arg(long)]
    clinvar: PathBuf,

    /// Path to PanelApp panel JSON (pre-downloaded).
    #[arg(long, conflicts_with = "panelapp_id")]
    panelapp: Option<PathBuf>,

    /// PanelApp panel ID to fetch live (e.g., 85 for hereditary neuropathy).
    #[arg(long, conflicts_with = "panelapp")]
    panelapp_id: Option<u32>,

    /// Output catalog JSON path.
    #[arg(long, default_value = "catalog.json")]
    output: PathBuf,

    /// Minimum ClinVar review stars (0-4). Default: 1 (criteria provided).
    #[arg(long, default_value = "1")]
    min_stars: u8,
}

/// Gene info extracted from PanelApp.
struct PanelGene {
    symbol: String,
    inheritance: Inheritance,
    confidence: u8, // 1=red, 2=amber, 3=green
}

fn main() {
    let cli = Cli::parse();

    // Step 1: Load PanelApp gene list
    let genes = if let Some(path) = &cli.panelapp {
        eprintln!("Loading PanelApp panel from {}", path.display());
        let data = std::fs::read_to_string(path).expect("failed to read PanelApp JSON");
        parse_panelapp_json(&data)
    } else if let Some(id) = cli.panelapp_id {
        eprintln!("Fetching PanelApp panel {id} from API...");
        let data = fetch_panelapp(id);
        parse_panelapp_json(&data)
    } else {
        eprintln!("ERROR: specify --panelapp or --panelapp-id");
        std::process::exit(1);
    };

    // Genes excluded from screening: require specialized callers
    // that cannot run in a browser WASM environment.
    const EXCLUDED_PARALOG: &[&str] = &["SMN1", "SMN2"];
    const EXCLUDED_REASON: &[(ExclusionReason, &str)] = &[
        (ExclusionReason::ParalogMappability,
         "Requires specialized caller (e.g. SMNCopyNumberCaller) due to near-identical paralog. Request dedicated SMA carrier screening if clinically indicated."),
        (ExclusionReason::MitochondrialHeteroplasmy,
         "Mitochondrial variants require heteroplasmy-aware calling. Request dedicated mtDNA panel if clinically indicated."),
    ];

    // Only include green + amber confidence genes
    let all_panel_genes: HashSet<String> = genes
        .iter()
        .filter(|g| g.confidence >= 2)
        .map(|g| g.symbol.clone())
        .collect();

    // Build exclusion list
    let mut excluded_genes: Vec<ExcludedGene> = Vec::new();
    let mut excluded_set: HashSet<String> = HashSet::new();

    for &gene in EXCLUDED_PARALOG {
        if all_panel_genes.contains(gene) {
            excluded_genes.push(ExcludedGene {
                symbol: gene.to_string(),
                reason: EXCLUDED_REASON[0].0,
                recommendation: EXCLUDED_REASON[0].1.to_string(),
            });
            excluded_set.insert(gene.to_string());
        }
    }
    // Exclude all mitochondrial genes (MT-* prefix)
    for gene in &all_panel_genes {
        if gene.starts_with("MT-") {
            excluded_genes.push(ExcludedGene {
                symbol: gene.clone(),
                reason: EXCLUDED_REASON[1].0,
                recommendation: EXCLUDED_REASON[1].1.to_string(),
            });
            excluded_set.insert(gene.clone());
        }
    }

    // Active gene set = panel genes minus excluded
    let gene_set: HashSet<String> = all_panel_genes.difference(&excluded_set).cloned().collect();
    let gene_inheritance: HashMap<String, Inheritance> = genes
        .into_iter()
        .filter(|g| g.confidence >= 2 && !excluded_set.contains(&g.symbol))
        .map(|g| (g.symbol.clone(), g.inheritance))
        .collect();

    eprintln!(
        "PanelApp: {} genes active, {} excluded",
        gene_set.len(),
        excluded_genes.len()
    );
    for eg in &excluded_genes {
        eprintln!("  excluded: {} ({:?})", eg.symbol, eg.reason);
    }

    // Step 2: Stream ClinVar VCF and filter
    eprintln!("Streaming ClinVar VCF: {}", cli.clinvar.display());
    let file = std::fs::File::open(&cli.clinvar).expect("failed to open ClinVar VCF");
    let reader: Box<dyn BufRead> = if cli.clinvar.extension().is_some_and(|e| e == "gz") {
        // ClinVar VCF is BGZF (blocked gzip). MultiGzDecoder handles
        // the concatenated gzip blocks correctly.
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut variants = Vec::new();
    let mut total_lines = 0u64;
    let mut skipped_gene = 0u64;
    let mut skipped_sig = 0u64;
    let mut skipped_stars = 0u64;

    for line in reader.lines() {
        let line = line.expect("failed to read VCF line");
        if line.starts_with('#') {
            continue;
        }
        total_lines += 1;
        if total_lines.is_multiple_of(500_000) {
            eprintln!(
                "  processed {} variants, {} kept...",
                total_lines,
                variants.len()
            );
        }

        if let Some(v) = parse_clinvar_line(&line, &gene_set, &gene_inheritance, cli.min_stars) {
            variants.push(v);
        } else {
            // Count skip reasons (approximate - just for logging)
            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() >= 8 {
                let info = fields[7];
                let gene = extract_info_field(info, "GENEINFO")
                    .and_then(|g| g.split(':').next().map(|s| s.to_string()));
                if let Some(g) = gene {
                    if !gene_set.contains(&g) {
                        skipped_gene += 1;
                        continue;
                    }
                }
                let clnsig = extract_info_field(info, "CLNSIG").unwrap_or_default();
                if !clnsig.contains("Pathogenic") {
                    skipped_sig += 1;
                    continue;
                }
                skipped_stars += 1;
            }
        }
    }

    // Sort by chrom + pos for efficient querying
    variants.sort_by(|a, b| a.chrom.cmp(&b.chrom).then_with(|| a.pos.cmp(&b.pos)));

    // Deduplicate by (chrom, pos, ref, alt)
    variants.dedup_by(|a, b| {
        a.chrom == b.chrom
            && a.pos == b.pos
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
    });

    eprintln!();
    eprintln!("ClinVar: {} total variants processed", total_lines);
    eprintln!("  skipped (gene not in panel): {}", skipped_gene);
    eprintln!("  skipped (not P/LP): {}", skipped_sig);
    eprintln!(
        "  skipped (below {} stars): {}",
        cli.min_stars, skipped_stars
    );
    eprintln!("  kept: {} pathogenic variants", variants.len());

    // Distinct genes
    let mut found_genes: Vec<&str> = variants.iter().map(|v| v.gene.as_str()).collect();
    found_genes.sort_unstable();
    found_genes.dedup();
    eprintln!("  across {} genes", found_genes.len());

    // Step 3: Build catalog
    let catalog = VariantCatalog {
        meta: CatalogMeta {
            generated_date: today_utc(),
            clinvar_date: clinvar_date_from_path(&cli.clinvar),
            panelapp_panel: format!(
                "PanelApp panel {} (Hereditary neuropathy)",
                cli.panelapp_id.unwrap_or(85)
            ),
            gene_count: found_genes.len(),
            variant_count: variants.len(),
            min_review_stars: cli.min_stars,
            excluded_genes,
        },
        variants,
    };

    // Step 4: Write output
    let json = serde_json::to_string_pretty(&catalog).expect("failed to serialize catalog");
    std::fs::write(&cli.output, &json).expect("failed to write output");
    eprintln!(
        "\nCatalog written to {} ({} bytes, {} variants across {} genes)",
        cli.output.display(),
        json.len(),
        catalog.meta.variant_count,
        catalog.meta.gene_count,
    );
}

fn parse_clinvar_line(
    line: &str,
    gene_set: &HashSet<String>,
    gene_inheritance: &HashMap<String, Inheritance>,
    min_stars: u8,
) -> Option<PathogenicVariant> {
    let fields: Vec<&str> = line.splitn(9, '\t').collect();
    if fields.len() < 8 {
        return None;
    }
    let chrom = fields[0];
    let pos: u32 = fields[1].parse().ok()?;
    let clinvar_id: u32 = fields[2].parse().ok()?;
    let ref_allele = fields[3];
    let alt_allele = fields[4];
    let info = fields[7];

    // Extract gene symbol from GENEINFO field (format: GENE:GENEID|GENE2:GENEID2)
    let gene = extract_info_field(info, "GENEINFO")
        .and_then(|g| g.split(':').next().map(|s| s.to_string()))?;

    if !gene_set.contains(&gene) {
        return None;
    }

    // Check clinical significance. ClinVar uses values like "Pathogenic",
    // "Likely_pathogenic", "Pathogenic/Likely_pathogenic", etc.
    // We accept any entry that contains "athogenic" (catches P, LP, P/LP).
    let clnsig = extract_info_field(info, "CLNSIG")?;
    if !clnsig.contains("athogenic") {
        return None;
    }
    // Exclude "Conflicting" entries that only have benign + pathogenic
    if clnsig.contains("Conflicting") {
        return None;
    }
    let significance = if clnsig == "Pathogenic" {
        ClinSig::Pathogenic
    } else {
        // LP, P/LP, and any other pathogenic-containing value
        ClinSig::LikelyPathogenic
    };

    // Check review status (CLNREVSTAT)
    let revstat = extract_info_field(info, "CLNREVSTAT").unwrap_or_default();
    let stars = review_stars(revstat);
    if stars < min_stars {
        return None;
    }

    // Extract condition
    let condition = extract_info_field(info, "CLNDN")
        .unwrap_or_default()
        .replace('_', " ");

    // Extract molecular consequence
    let consequence = extract_info_field(info, "MC")
        .and_then(|mc| mc.split('|').nth(1).map(|s| s.to_string()))
        .unwrap_or_default();

    let inheritance = gene_inheritance
        .get(&gene)
        .copied()
        .unwrap_or(Inheritance::Unknown);

    Some(PathogenicVariant {
        chrom: chrom.to_string(),
        pos,
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
        gene,
        clinvar_id,
        significance,
        review_stars: stars,
        condition,
        inheritance,
        consequence,
    })
}

fn extract_info_field<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    let prefix = format!("{key}=");
    info.split(';')
        .find(|f| f.starts_with(&prefix))
        .map(|f| &f[prefix.len()..])
}

fn review_stars(revstat: &str) -> u8 {
    match revstat {
        s if s.contains("practice_guideline") => 4,
        s if s.contains("reviewed_by_expert_panel") => 3,
        s if s.contains("criteria_provided,_multiple_submitters,_no_conflicts") => 2,
        s if s.contains("criteria_provided,_conflicting") => 1,
        s if s.contains("criteria_provided,_single_submitter") => 1,
        _ => 0,
    }
}

fn parse_panelapp_json(json: &str) -> Vec<PanelGene> {
    let data: serde_json::Value = serde_json::from_str(json).expect("invalid PanelApp JSON");
    let genes = data["genes"]
        .as_array()
        .expect("no genes array in PanelApp JSON");

    genes
        .iter()
        .filter_map(|g| {
            let symbol = g["gene_data"]["gene_symbol"].as_str()?.to_string();
            let confidence: u8 = g["confidence_level"]
                .as_str()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
            let moi = g["mode_of_inheritance"].as_str().unwrap_or("");
            let inheritance = parse_inheritance(moi);
            Some(PanelGene {
                symbol,
                inheritance,
                confidence,
            })
        })
        .collect()
}

fn parse_inheritance(moi: &str) -> Inheritance {
    let lower = moi.to_lowercase();
    let has_ad = lower.contains("monoallelic");
    let has_ar = lower.contains("biallelic");
    let has_x = lower.contains("x-linked");

    if has_x {
        Inheritance::XLinked
    } else if has_ad && has_ar {
        Inheritance::Both
    } else if has_ar {
        Inheritance::AutosomalRecessive
    } else if has_ad {
        Inheritance::AutosomalDominant
    } else {
        Inheritance::Unknown
    }
}

fn fetch_panelapp(panel_id: u32) -> String {
    let url =
        format!("https://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/?format=json");
    eprintln!("  GET {url}");

    let output = std::process::Command::new("curl")
        .args(["-s", "-f", &url])
        .output()
        .expect("failed to run curl");

    if !output.status.success() {
        eprintln!(
            "ERROR: PanelApp API returned status {}",
            output.status.code().unwrap_or(-1)
        );
        std::process::exit(1);
    }

    String::from_utf8(output.stdout).expect("PanelApp response is not UTF-8")
}

fn today_utc() -> String {
    format_ymd(std::time::SystemTime::now())
}

fn clinvar_date_from_path(path: &PathBuf) -> String {
    std::fs::metadata(path)
        .and_then(|m| m.modified())
        .map(format_ymd)
        .unwrap_or_else(|_| "unknown".to_string())
}

/// Format a `SystemTime` as `YYYY-MM-DD` in UTC using a proper
/// proleptic-Gregorian calculation. Accurate across leap years.
fn format_ymd(t: std::time::SystemTime) -> String {
    let secs = t
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let days = secs.div_euclid(86_400);
    let (y, m, d) = days_to_ymd(days);
    format!("{y:04}-{m:02}-{d:02}")
}

/// Convert days since 1970-01-01 (UTC) to `(year, month, day)`.
fn days_to_ymd(days_since_epoch: i64) -> (i32, u32, u32) {
    // Shift to days since 0000-03-01 so the "leap-year exceptions" align
    // with the start of a 400-year cycle. The algorithm is from Howard
    // Hinnant's "chrono-Compatible Low-Level Date Algorithms".
    let z = days_since_epoch + 719_468;
    let era = z.div_euclid(146_097);
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = (yoe as i64) + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y as i32, m, d)
}

#[cfg(test)]
mod tests {
    use super::days_to_ymd;

    #[test]
    fn days_to_ymd_known_dates() {
        assert_eq!(days_to_ymd(0), (1970, 1, 1));
        assert_eq!(days_to_ymd(59), (1970, 3, 1));
        assert_eq!(days_to_ymd(365), (1971, 1, 1));
        // 2000 was a leap year; Feb 29 was day 11_017.
        assert_eq!(days_to_ymd(11_016), (2000, 2, 29));
        assert_eq!(days_to_ymd(11_017), (2000, 3, 1));
        // 2026-04-20
        assert_eq!(days_to_ymd(20_563), (2026, 4, 20));
    }
}
