//! WASM-side pileup implementation for BAM and CRAM.

use nds_core::catalog::{PathogenicVariant, VariantCatalog};
use nds_core::pileup::{build_pileup_result, CheckedVariant, PileupConfig, PileupResult};
use noodles_sam::alignment::record::cigar::op::Kind;
use std::collections::HashMap;
use std::io::{BufReader, Cursor};

// ---------------------------------------------------------------------------
// BAM base extraction (lazy record API - concrete types, no trait objects)
// ---------------------------------------------------------------------------

fn base_at_pos_bam(record: &noodles_bam::Record, target_pos: u32, min_bq: u8) -> Option<u8> {
    let aln_start = record.alignment_start()?.ok()?.get() as u32;
    let cigar = record.cigar();
    let sequence = record.sequence();
    let quals = record.quality_scores();

    let mut ref_pos = aln_start;
    let mut read_pos: usize = 0;

    for op_result in cigar.iter() {
        let op = op_result.ok()?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for _ in 0..len {
                    if ref_pos == target_pos {
                        let base = sequence.get(read_pos)?;
                        // Check base quality (BAM lazy QualityScores)
                        let bq = quals.iter().nth(read_pos).unwrap_or(0);
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

// ---------------------------------------------------------------------------
// Classify variants from counts
// ---------------------------------------------------------------------------

/// Key for per-variant counting: position + specific alt allele.
type VariantKey = (u32, u8); // (pos, alt_base)

/// Per-variant alt allele counts with strand tracking.
#[derive(Default, Clone)]
struct AlleleCount {
    alt: u32,
    forward_alt: u32,
    reverse_alt: u32,
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

// ---------------------------------------------------------------------------
// BAM pileup
// ---------------------------------------------------------------------------

pub fn pileup_bam(
    bam_bytes: &[u8],
    bai_bytes: &[u8],
    catalog: &VariantCatalog,
    config: &PileupConfig,
) -> Result<PileupResult, String> {
    let index = noodles_bam::bai::io::Reader::new(Cursor::new(bai_bytes))
        .read_index()
        .map_err(|e| format!("BAI parse error: {e}"))?;
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_reader(BufReader::new(Cursor::new(bam_bytes)))
        .map_err(|e| format!("BAM reader build error: {e}"))?;
    let header = reader
        .read_header()
        .map_err(|e| format!("BAM header error: {e}"))?;
    run_pileup_bam_reader(&mut reader, &header, catalog, config)
}

// ---------------------------------------------------------------------------
// CRAM pileup
// ---------------------------------------------------------------------------

pub fn pileup_cram(
    cram_bytes: &[u8],
    crai_bytes: &[u8],
    ref_seqs: Vec<(String, Vec<u8>)>,
    catalog: &VariantCatalog,
    config: &PileupConfig,
    sex: nds_core::Sex,
) -> Result<PileupResult, String> {
    let crai_index = noodles_cram::crai::io::Reader::new(Cursor::new(crai_bytes))
        .read_index()
        .map_err(|e| format!("CRAI parse error: {e}"))?;

    let mut reader = build_cram_reader(Cursor::new(cram_bytes), ref_seqs);
    let header = reader
        .read_header()
        .map_err(|e| format!("CRAM header error: {e}"))?;

    nds_core::pileup_cram::pileup_cram(&mut reader, &header, &crai_index, catalog, config, sex)
}

fn build_cram_reader<R: std::io::Read>(
    inner: R,
    ref_seqs: Vec<(String, Vec<u8>)>,
) -> noodles_cram::io::Reader<BufReader<R>> {
    let mut builder = noodles_cram::io::reader::Builder::default();
    if !ref_seqs.is_empty() {
        let records: Vec<noodles_fasta::Record> = ref_seqs
            .into_iter()
            .map(|(name, seq_bytes)| {
                let def = noodles_fasta::record::Definition::new(name, None);
                let seq = noodles_fasta::record::Sequence::from(seq_bytes);
                noodles_fasta::Record::new(def, seq)
            })
            .collect();
        let repo = noodles_fasta::repository::Repository::new(records);
        builder = builder.set_reference_sequence_repository(repo);
    }
    builder.build_from_reader(BufReader::new(inner))
}

// ---------------------------------------------------------------------------
// Sparse BAM pileup (from File.slice() data via OffsetCursor)
// ---------------------------------------------------------------------------

pub fn pileup_bam_sparse(
    data: Vec<u8>,
    range_descs: &[(u64, u64)],
    bai_bytes: Vec<u8>,
    catalog: &VariantCatalog,
    config: &PileupConfig,
) -> Result<PileupResult, String> {
    use nds_core::sources::offset_cursor::OffsetCursor;

    let cursor = OffsetCursor::from_concatenated(data, range_descs);
    let index = noodles_bam::bai::io::Reader::new(Cursor::new(bai_bytes))
        .read_index()
        .map_err(|e| format!("BAI parse error: {e}"))?;

    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_reader(BufReader::new(cursor))
        .map_err(|e| format!("BAM reader build error: {e}"))?;
    let header = reader
        .read_header()
        .map_err(|e| format!("BAM header error: {e}"))?;

    run_pileup_bam_reader(&mut reader, &header, catalog, config)
}

// ---------------------------------------------------------------------------
// Sparse CRAM pileup (from File.slice() data via OffsetCursor + reference)
// ---------------------------------------------------------------------------

pub fn pileup_cram_sparse(
    data: Vec<u8>,
    range_descs: &[(u64, u64)],
    crai_bytes: Vec<u8>,
    ref_seqs: Vec<(String, Vec<u8>)>,
    catalog: &VariantCatalog,
    config: &PileupConfig,
    sex: nds_core::Sex,
) -> Result<PileupResult, String> {
    use nds_core::sources::offset_cursor::OffsetCursor;

    let crai_index = noodles_cram::crai::io::Reader::new(Cursor::new(crai_bytes))
        .read_index()
        .map_err(|e| format!("CRAI parse error: {e}"))?;

    let cursor = OffsetCursor::from_concatenated(data, range_descs);
    let mut reader = build_cram_reader(cursor, ref_seqs);
    let header = reader
        .read_header()
        .map_err(|e| format!("CRAM header error: {e}"))?;

    nds_core::pileup_cram::pileup_cram(&mut reader, &header, &crai_index, catalog, config, sex)
}

/// Shared BAM pileup logic over any indexed BAM reader.
fn run_pileup_bam_reader<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles_bam::io::IndexedReader<noodles_bgzf::io::Reader<R>>,
    header: &noodles_sam::Header,
    catalog: &VariantCatalog,
    config: &PileupConfig,
) -> Result<PileupResult, String> {
    let mut by_chrom: HashMap<&str, Vec<&PathogenicVariant>> = HashMap::new();
    for v in &catalog.variants {
        if v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
            by_chrom.entry(v.chrom.as_str()).or_default().push(v);
        }
    }

    let mut all_checked: Vec<CheckedVariant> = Vec::new();
    // Add non-SNV variants as NoCoverage (indel pileup not yet supported)
    for v in &catalog.variants {
        if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
            all_checked.push(nds_core::pileup_cram::no_coverage_variant(v));
        }
    }

    for (catalog_chrom, variants) in &by_chrom {
        let header_chrom = match nds_core::pileup_cram::resolve_chrom(catalog_chrom, header) {
            Some(c) => c,
            None => {
                for v in variants {
                    all_checked.push(nds_core::pileup_cram::no_coverage_variant(v));
                }
                continue;
            }
        };

        let mut alt_counts: HashMap<VariantKey, AlleleCount> = HashMap::new();
        let mut ref_counts: HashMap<u32, u32> = HashMap::new();
        let mut ref_bases: HashMap<u32, u8> = HashMap::new();
        let mut variant_positions: Vec<u32> = Vec::new();
        for v in variants {
            if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
                continue;
            }
            alt_counts.insert((v.pos, v.alt_allele.as_bytes()[0]), AlleleCount::default());
            ref_bases.insert(v.pos, v.ref_allele.as_bytes()[0]);
            variant_positions.push(v.pos);
        }
        variant_positions.sort_unstable();
        variant_positions.dedup();
        if variant_positions.is_empty() {
            all_checked.extend(
                variants
                    .iter()
                    .map(|v| nds_core::pileup_cram::no_coverage_variant(v)),
            );
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
            .query(header, &region)
            .map_err(|e| format!("BAM query error on {header_chrom}: {e}"))?;
        for record_result in query.records() {
            let record = match record_result {
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
                Some(Ok(pos)) => pos.get() as u32,
                _ => continue,
            };
            let is_reverse = flags.is_reverse_complemented();
            for &pos in &variant_positions {
                if pos >= aln_start {
                    let base = base_at_pos_bam(&record, pos, config.min_base_qual);
                    count_base_hit(&mut alt_counts, pos, base, is_reverse);
                    if let Some(&rb) = ref_bases.get(&pos) {
                        count_ref_at_pos(&mut ref_counts, pos, base, rb);
                    }
                }
            }
        }
        for v in variants {
            if v.ref_allele.len() != 1 || v.alt_allele.len() != 1 {
                all_checked.push(nds_core::pileup_cram::no_coverage_variant(v));
                continue;
            }
            let alt_base = v.alt_allele.as_bytes()[0];
            let ac = alt_counts
                .get(&(v.pos, alt_base))
                .cloned()
                .unwrap_or_default();
            let (alt_count, fwd_alt, rev_alt) = (ac.alt, ac.forward_alt, ac.reverse_alt);
            let ref_count = ref_counts.get(&v.pos).cloned().unwrap_or(0);
            let call = nds_core::pileup::classify_pileup_with_strand(
                ref_count, alt_count, fwd_alt, rev_alt, config,
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
