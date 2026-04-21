//! CRAM-based pileup genotyping.
//!
//! Runs the CRAM → pileup → classification pipeline on any
//! `Read + Seek` source, so both the WASM entry point and native tests
//! exercise the same code.

use std::collections::HashMap;
use std::io::{Read, Seek};

use noodles_sam::alignment::record::cigar::op::Kind;

use crate::catalog::{PathogenicVariant, VariantCatalog};
use crate::pileup::{
    build_pileup_result, classify_pileup_with_context, CheckedVariant, PileupConfig, PileupResult,
    VariantCall,
};
use crate::Sex;

/// Key for per-variant counting: (position, specific alt base).
type VariantKey = (u32, u8);

/// Per-variant alt allele counts with strand tracking.
#[derive(Default, Clone)]
struct AlleleCount {
    alt: u32,
    forward_alt: u32,
    reverse_alt: u32,
}

/// Extract base at a reference position from a CRAM RecordBuf.
/// Returns None if position is deleted, outside read, or base quality < min_bq.
fn base_at_pos_cram(
    record: &noodles_sam::alignment::RecordBuf,
    target_pos: u32,
    min_bq: u8,
) -> Option<u8> {
    let aln_start = record.alignment_start()?.get() as u32;
    let cigar = record.cigar();
    let sequence = record.sequence();
    let quals = record.quality_scores();

    let mut ref_pos = aln_start;
    let mut read_pos: usize = 0;

    for op in cigar.as_ref() {
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if ref_pos == target_pos {
                        let base = sequence.as_ref().get(read_pos).cloned()?;
                        let bq = quals.as_ref().get(read_pos).cloned().unwrap_or(0);
                        if bq < min_bq {
                            return None;
                        }
                        return Some(base.to_ascii_uppercase());
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Kind::Insertion | Kind::SoftClip => read_pos += len,
            Kind::Deletion | Kind::Skip => {
                if target_pos >= ref_pos && target_pos < ref_pos + len as u32 {
                    return None;
                }
                ref_pos += len as u32;
            }
            Kind::HardClip | Kind::Pad => {}
        }
        if ref_pos > target_pos {
            break;
        }
    }
    None
}

fn count_base_hit(
    counts: &mut HashMap<VariantKey, AlleleCount>,
    pos: u32,
    base: Option<u8>,
    is_reverse: bool,
) {
    if let Some(base) = base {
        if let Some(entry) = counts.get_mut(&(pos, base)) {
            entry.alt += 1;
            if is_reverse {
                entry.reverse_alt += 1;
            } else {
                entry.forward_alt += 1;
            }
        }
    }
}

fn count_ref_at_pos(ref_counts: &mut HashMap<u32, u32>, pos: u32, base: Option<u8>, ref_base: u8) {
    if let Some(base) = base {
        if base == ref_base {
            *ref_counts.entry(pos).or_default() += 1;
        }
    }
}

/// Resolve a ClinVar chromosome name (e.g. "1", "X") to the name used
/// in the BAM/CRAM header (e.g. "chr1", "chrX").
pub fn resolve_chrom(chrom: &str, header: &noodles_sam::Header) -> Option<String> {
    let refs = header.reference_sequences();
    if refs.contains_key(chrom.as_bytes()) {
        return Some(chrom.to_string());
    }
    let with_chr = format!("chr{chrom}");
    if refs.contains_key(with_chr.as_bytes()) {
        return Some(with_chr);
    }
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if refs.contains_key(stripped.as_bytes()) {
            return Some(stripped.to_string());
        }
    }
    None
}

pub fn no_coverage_variant(v: &PathogenicVariant) -> CheckedVariant {
    CheckedVariant {
        gene: v.gene.clone(),
        chrom: v.chrom.clone(),
        pos: v.pos,
        ref_allele: v.ref_allele.clone(),
        alt_allele: v.alt_allele.clone(),
        clinvar_id: v.clinvar_id,
        significance: v.significance,
        review_stars: v.review_stars,
        condition: v.condition.clone(),
        inheritance: v.inheritance,
        consequence: v.consequence.clone(),
        call: VariantCall::NoCoverage,
    }
}

/// Run pileup genotyping on a CRAM reader against a variant catalog.
///
/// This is the core pileup pipeline: for each chromosome with catalog
/// variants, query the CRAM, extract bases at variant positions, count
/// alt/ref with strand tracking, and classify.
pub fn pileup_cram<R: Read + Seek>(
    reader: &mut noodles_cram::io::Reader<R>,
    header: &noodles_sam::Header,
    crai_index: &noodles_cram::crai::Index,
    catalog: &VariantCatalog,
    config: &PileupConfig,
    sex: Sex,
) -> Result<PileupResult, String> {
    let mut by_chrom: HashMap<&str, Vec<&PathogenicVariant>> = HashMap::new();
    for v in &catalog.variants {
        if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
            by_chrom.entry(v.chrom.as_str()).or_default().push(v);
        }
    }

    let mut all_checked: Vec<CheckedVariant> = Vec::new();
    for v in &catalog.variants {
        if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
            all_checked.push(no_coverage_variant(v));
        }
    }

    for (catalog_chrom, variants) in &by_chrom {
        let header_chrom = match resolve_chrom(catalog_chrom, header) {
            Some(c) => c,
            None => {
                for v in variants {
                    all_checked.push(no_coverage_variant(v));
                }
                continue;
            }
        };

        let mut alt_counts: HashMap<VariantKey, AlleleCount> = HashMap::new();
        let mut ref_counts: HashMap<u32, u32> = HashMap::new();
        let mut ref_bases: HashMap<u32, u8> = HashMap::new();
        let mut variant_positions: Vec<u32> = Vec::new();

        for v in variants {
            alt_counts.insert((v.pos, v.alt_allele.as_bytes()[0]), AlleleCount::default());
            ref_bases.insert(v.pos, v.ref_allele.as_bytes()[0]);
            variant_positions.push(v.pos);
        }
        variant_positions.sort_unstable();
        variant_positions.dedup();

        if variant_positions.is_empty() {
            all_checked.extend(variants.iter().map(|v| no_coverage_variant(v)));
            continue;
        }

        let min_pos = *variant_positions.first().unwrap();
        let max_pos = *variant_positions.last().unwrap();
        let start = noodles_core::Position::try_from(min_pos as usize)
            .map_err(|e| format!("invalid start: {e}"))?;
        let end = noodles_core::Position::try_from(max_pos as usize)
            .map_err(|e| format!("invalid end: {e}"))?;
        let region = noodles_core::Region::new(header_chrom.clone(), start..=end);

        let query = reader
            .query(header, crai_index, &region)
            .map_err(|e| format!("CRAM query error on {header_chrom}: {e}"))?;

        for record_result in query {
            let record: noodles_sam::alignment::RecordBuf = match record_result {
                Ok(r) => r,
                Err(_) => continue,
            };
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let mapq: u8 = record.mapping_quality().map(u8::from).unwrap_or(0);
            if mapq < config.min_mapq {
                continue;
            }
            let aln_start = match record.alignment_start() {
                Some(pos) => pos.get() as u32,
                None => continue,
            };

            let ref_span: u32 = record
                .cigar()
                .as_ref()
                .iter()
                .map(|op| match op.kind() {
                    Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion
                    | Kind::Skip => op.len() as u32,
                    _ => 0,
                })
                .sum();
            let aln_end = aln_start + ref_span.saturating_sub(1);

            let is_reverse = flags.is_reverse_complemented();
            for &pos in &variant_positions {
                if pos >= aln_start && pos <= aln_end {
                    let base = base_at_pos_cram(&record, pos, config.min_base_qual);
                    count_base_hit(&mut alt_counts, pos, base, is_reverse);
                    if let Some(&rb) = ref_bases.get(&pos) {
                        count_ref_at_pos(&mut ref_counts, pos, base, rb);
                    }
                }
            }
        }

        for v in variants {
            let alt_base = v.alt_allele.as_bytes()[0];
            let ac = alt_counts
                .get(&(v.pos, alt_base))
                .cloned()
                .unwrap_or_default();
            let (alt_count, fwd_alt, rev_alt) = (ac.alt, ac.forward_alt, ac.reverse_alt);
            let ref_count = ref_counts.get(&v.pos).cloned().unwrap_or(0);

            let call = classify_pileup_with_context(
                ref_count,
                alt_count,
                fwd_alt,
                rev_alt,
                config,
                &header_chrom,
                sex,
            );
            all_checked.push(CheckedVariant {
                gene: v.gene.clone(),
                chrom: v.chrom.clone(),
                pos: v.pos,
                ref_allele: v.ref_allele.clone(),
                alt_allele: v.alt_allele.clone(),
                clinvar_id: v.clinvar_id,
                significance: v.significance,
                review_stars: v.review_stars,
                condition: v.condition.clone(),
                inheritance: v.inheritance,
                consequence: v.consequence.clone(),
                call,
            });
        }
    }

    Ok(build_pileup_result(all_checked, catalog))
}
