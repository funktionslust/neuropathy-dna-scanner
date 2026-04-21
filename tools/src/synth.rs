//! Synthetic BAM generation engine.
//!
//! Produces deterministic, coordinate-sorted BAM files at specified
//! coverage profiles. Used by the `synth_bam` binary and by the xtask
//! fixture generator.

use std::io;
use std::num::NonZeroUsize;
use std::path::Path;

use noodles_sam::header::record::value::map::reference_sequence::tag;
use noodles_sam::header::record::value::{map::ReferenceSequence, Map};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

/// Read length for all synthetic reads.
const READ_LENGTH: u32 = 150;

/// GRCh38 standard chromosome names and lengths (chr1-chr22, chrX, chrY).
const GRCH38_CHROMS: &[(&str, u64)] = &[
    ("chr1", 248_956_422),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
    ("chr6", 170_805_979),
    ("chr7", 159_345_973),
    ("chr8", 145_138_636),
    ("chr9", 138_394_717),
    ("chr10", 133_797_422),
    ("chr11", 135_086_622),
    ("chr12", 133_275_309),
    ("chr13", 114_364_328),
    ("chr14", 107_043_718),
    ("chr15", 101_991_189),
    ("chr16", 90_338_345),
    ("chr17", nds_core::GRCH38_CHR17_LENGTH),
    ("chr18", 80_373_285),
    ("chr19", 58_617_616),
    ("chr20", 64_444_167),
    ("chr21", 46_709_983),
    ("chr22", 50_818_468),
    ("chrX", 156_040_895),
    ("chrY", 57_227_415),
];

/// The real chr17 M5 tag from GRCh38 (validated against real WGS CRAMs).
const CHR17_M5: &str = "f9a0fb01553adb183568e3eb9d8626db";

/// A parsed coverage specification: chromosome, start, end, depth.
#[derive(Debug, Clone)]
pub struct CoverageSpec {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub depth: u32,
}

impl CoverageSpec {
    /// Parse from the CLI format `chr:start-end=depth`.
    pub fn parse(s: &str) -> Result<Self, String> {
        let (region, depth_str) = s
            .rsplit_once('=')
            .ok_or_else(|| format!("missing '=' in coverage spec: {s}"))?;
        let depth: u32 = depth_str
            .parse()
            .map_err(|_| format!("invalid depth: {depth_str}"))?;
        let (chrom, range) = region
            .split_once(':')
            .ok_or_else(|| format!("missing ':' in region: {region}"))?;
        let (start_str, end_str) = range
            .split_once('-')
            .ok_or_else(|| format!("missing '-' in range: {range}"))?;
        let start: u32 = start_str
            .parse()
            .map_err(|_| format!("invalid start: {start_str}"))?;
        let end: u32 = end_str
            .parse()
            .map_err(|_| format!("invalid end: {end_str}"))?;
        if start == 0 || end == 0 || start > end {
            return Err(format!("invalid range: {start}-{end}"));
        }
        Ok(Self {
            chrom: chrom.to_string(),
            start,
            end,
            depth,
        })
    }
}

/// Build a GRCh38 SAM header with 24 standard chromosomes.
pub fn build_grch38_header() -> noodles_sam::Header {
    let mut builder = noodles_sam::Header::builder();
    for &(name, length) in GRCH38_CHROMS {
        let len = NonZeroUsize::try_from(length as usize).unwrap();
        let mut ref_seq = Map::<ReferenceSequence>::new(len);
        // Add M5 tag: real value for chr17, deterministic placeholder for others.
        let m5 = if name == "chr17" {
            CHR17_M5.to_string()
        } else {
            format!("{:032x}", md5_placeholder(name))
        };
        ref_seq
            .other_fields_mut()
            .insert(tag::MD5_CHECKSUM, m5.into());
        builder = builder.add_reference_sequence(name, ref_seq);
    }
    builder.build()
}

/// Deterministic 128-bit placeholder for non-chr17 M5 tags. Uses a simple
/// FNV-1a-style hash that is stable across Rust versions and platforms
/// (unlike `DefaultHasher`, which may change its algorithm between releases).
fn md5_placeholder(name: &str) -> u128 {
    // FNV-1a 64-bit on the name bytes, then again on the reversed bytes,
    // concatenated into 128 bits. This is NOT cryptographic  - it only needs
    // to be deterministic and collision-resistant enough for 24 chromosomes.
    fn fnv1a(data: &[u8]) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for &b in data {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        h
    }
    let a = fnv1a(name.as_bytes());
    let rev: Vec<u8> = name.bytes().rev().collect();
    let b = fnv1a(&rev);
    ((a as u128) << 64) | (b as u128)
}

/// A variant to spike into the synthetic reads. At the specified position,
/// a fraction of overlapping reads will carry the alt allele instead of
/// random bases. Used to create positive-control fixtures for SNV testing.
/// Strand override for spike reads.
#[derive(Debug, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone)]
pub struct VariantSpike {
    pub chrom: String,
    /// 1-based position on the reference.
    pub pos: u32,
    /// Reference allele (single base for SNVs).
    pub ref_base: u8,
    /// Alt allele to spike in.
    pub alt_base: u8,
    /// Fraction of reads that should carry the alt allele (0.0-1.0).
    /// 0.5 = heterozygous, 1.0 = homozygous alt.
    pub allele_fraction: f64,
    /// Force all alt-carrying reads to one strand (for strand-bias testing).
    pub strand_override: Option<Strand>,
    /// Override base quality at the spike position for alt reads.
    pub alt_base_quality: Option<u8>,
    /// Override MAPQ for reads carrying the alt allele.
    pub alt_mapq: Option<u8>,
}

impl VariantSpike {
    /// Create a simple het spike with no quality/strand overrides.
    pub fn het(chrom: impl Into<String>, pos: u32, ref_base: u8, alt_base: u8) -> Self {
        Self {
            chrom: chrom.into(),
            pos,
            ref_base,
            alt_base,
            allele_fraction: 0.5,
            strand_override: None,
            alt_base_quality: None,
            alt_mapq: None,
        }
    }
}

/// A single synthetic read to be written to the BAM/CRAM.
pub struct SynthRead {
    pub ref_id: usize,
    pub position: u32, // 1-based
    pub bases: Vec<u8>,
    pub quality_scores: Vec<u8>,
    pub mapq: u8,
    pub is_reverse: bool,
}

/// Generate reads for a set of coverage specs with optional variant spikes.
/// Reads overlapping spike positions have their base at that position set
/// to the alt allele with probability = allele_fraction, or the ref allele
/// otherwise (instead of a random base).
pub fn generate_reads(
    header: &noodles_sam::Header,
    specs: &[CoverageSpec],
    seed: u64,
) -> Vec<SynthRead> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut reads = Vec::new();

    for spec in specs {
        let ref_id = header
            .reference_sequences()
            .get_index_of(spec.chrom.as_bytes())
            .unwrap_or_else(|| panic!("chromosome {} not in header", spec.chrom));

        let region_length = spec.end - spec.start + 1;
        let num_reads = (u64::from(spec.depth) * u64::from(region_length)) / u64::from(READ_LENGTH);

        for i in 0..num_reads {
            // Evenly spaced with small jitter.
            let base_pos = spec.start as u64 + (i * u64::from(region_length)) / num_reads;
            let jitter: i32 = rng.random_range(-50..=50);
            // Clamp to [spec.start, spec.end - READ_LENGTH + 1] so the
            // 150 bp read never extends past the region's end coordinate.
            let max_pos = spec.end.saturating_sub(READ_LENGTH - 1).max(spec.start);
            let position = (base_pos as i64 + i64::from(jitter))
                .max(spec.start as i64)
                .min(max_pos as i64) as u32;

            let bases: Vec<u8> = (0..READ_LENGTH)
                .map(|_| b"ACGT"[rng.random_range(0..4usize)])
                .collect();

            reads.push(SynthRead {
                ref_id,
                position,
                bases,
                quality_scores: vec![30; READ_LENGTH as usize],
                mapq: 60,
                is_reverse: rng.random_range(0..2u8) == 1,
            });
        }
    }

    reads.sort_by_key(|r| (r.ref_id, r.position));
    reads
}

/// Write a complete BAM file from the header and sorted reads.
pub fn write_bam(path: &Path, header: &noodles_sam::Header, reads: &[SynthRead]) -> io::Result<()> {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::alignment::RecordBuf;

    let file = std::fs::File::create(path)?;
    let mut writer = noodles_bam::io::Writer::new(file);
    writer.write_header(header)?;

    let cigar = Cigar::from(vec![Op::new(Kind::Match, READ_LENGTH as usize)]);
    let base_flags = Flags::SEGMENTED | Flags::MATE_UNMAPPED | Flags::FIRST_SEGMENT;

    for read in reads {
        let pos =
            noodles_core::Position::try_from(read.position as usize).map_err(io::Error::other)?;
        let mut flags = base_flags;
        if read.is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        let record = RecordBuf::builder()
            .set_name("synth")
            .set_reference_sequence_id(read.ref_id)
            .set_alignment_start(pos)
            .set_cigar(cigar.clone())
            .set_sequence(Sequence::from(read.bases.clone()))
            .set_quality_scores(QualityScores::from(read.quality_scores.clone()))
            .set_mapping_quality(MappingQuality::new(read.mapq).unwrap())
            .set_flags(flags)
            .build();
        writer.write_alignment_record(header, &record)?;
    }

    writer.try_finish()?;
    Ok(())
}

/// Write a CRAM file from header and sorted reads, with reference.
/// The reference is needed for CRAM's reference-based compression.
pub fn write_cram(
    path: &Path,
    header: &noodles_sam::Header,
    reads: &[SynthRead],
    reference_path: &Path,
) -> io::Result<()> {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::alignment::RecordBuf;

    // Build reference repository from indexed FASTA (.fai required)
    let fasta_reader = noodles_fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| io::Error::other(format!("open indexed FASTA: {e} (needs .fai)")))?;
    let repo_adapter = noodles_fasta::repository::adapters::IndexedReader::new(fasta_reader);
    let repo = noodles_fasta::Repository::new(repo_adapter);

    let file = std::fs::File::create(path)?;
    let mut writer = noodles_cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repo)
        .build_from_writer(file);
    writer.write_header(header)?;

    let cigar = Cigar::from(vec![Op::new(Kind::Match, READ_LENGTH as usize)]);
    let base_flags = Flags::SEGMENTED | Flags::MATE_UNMAPPED | Flags::FIRST_SEGMENT;

    for read in reads {
        let pos =
            noodles_core::Position::try_from(read.position as usize).map_err(io::Error::other)?;
        let mut flags = base_flags;
        if read.is_reverse {
            flags |= Flags::REVERSE_COMPLEMENTED;
        }
        let record = RecordBuf::builder()
            .set_name("synth")
            .set_reference_sequence_id(read.ref_id)
            .set_alignment_start(pos)
            .set_cigar(cigar.clone())
            .set_sequence(Sequence::from(read.bases.clone()))
            .set_quality_scores(QualityScores::from(read.quality_scores.clone()))
            .set_mapping_quality(MappingQuality::new(read.mapq).unwrap())
            .set_flags(flags)
            .build();
        writer.write_alignment_record(header, &record)?;
    }

    writer.try_finish(header)?;
    Ok(())
}

/// Indexed FASTA reader cached across calls for performance.
/// Without this, each read_ref_bases call reopens and re-indexes the 3 GB file.
use std::cell::RefCell;
thread_local! {
    static FASTA_CACHE: RefCell<Option<(std::path::PathBuf, noodles_fasta::io::IndexedReader<noodles_fasta::io::BufReader<std::fs::File>>)>> = const { RefCell::new(None) };
}

/// Read reference bases from an indexed FASTA (.fai required).
/// Returns READ_LENGTH bases starting at the 1-based position.
/// Uses a thread-local cache to avoid reopening the file per read.
pub fn read_ref_bases(fasta_path: &Path, chrom: &str, pos_1based: u32) -> io::Result<Vec<u8>> {
    FASTA_CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        // Open and cache the indexed reader on first call
        if cache.as_ref().is_none_or(|(p, _)| p != fasta_path) {
            let reader = noodles_fasta::io::indexed_reader::Builder::default()
                .build_from_path(fasta_path)
                .map_err(|e| io::Error::other(format!("open indexed FASTA: {e}")))?;
            *cache = Some((fasta_path.to_path_buf(), reader));
        }
        let (_, reader) = cache.as_mut().unwrap();

        let start = noodles_core::Position::try_from(pos_1based as usize)
            .map_err(|e| io::Error::other(format!("invalid pos: {e}")))?;
        let end = noodles_core::Position::try_from((pos_1based + READ_LENGTH - 1) as usize)
            .map_err(|e| io::Error::other(format!("invalid end: {e}")))?;
        let region = noodles_core::Region::new(chrom, start..=end);

        let record = reader
            .query(&region)
            .map_err(|e| io::Error::other(format!("FASTA query {chrom}:{pos_1based}: {e}")))?;

        let seq = record.sequence();
        let bytes: Vec<u8> = seq.as_ref().to_vec();

        if bytes.len() >= READ_LENGTH as usize {
            Ok(bytes[..READ_LENGTH as usize].to_vec())
        } else {
            let mut padded = bytes;
            padded.resize(READ_LENGTH as usize, b'N');
            Ok(padded)
        }
    })
}

/// Generate reads using actual reference bases (CRAM-friendly).
/// Each read's bases match the reference at its alignment position,
/// except at spike positions where the alt allele is injected.
pub fn generate_reads_with_ref(
    header: &noodles_sam::Header,
    specs: &[CoverageSpec],
    spikes: &[VariantSpike],
    fasta_path: &Path,
    seed: u64,
) -> Vec<SynthRead> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut reads = Vec::new();

    for spec in specs {
        let ref_id = header
            .reference_sequences()
            .get_index_of(spec.chrom.as_bytes())
            .unwrap_or_else(|| panic!("chromosome {} not in header", spec.chrom));

        let region_length = spec.end - spec.start + 1;
        let num_reads = (u64::from(spec.depth) * u64::from(region_length)) / u64::from(READ_LENGTH);

        for i in 0..num_reads {
            let base_pos = spec.start as u64 + (i * u64::from(region_length)) / num_reads;
            let jitter: i32 = rng.random_range(-50..=50);
            let max_pos = spec.end.saturating_sub(READ_LENGTH - 1).max(spec.start);
            let position = (base_pos as i64 + i64::from(jitter))
                .max(spec.start as i64)
                .min(max_pos as i64) as u32;

            // Read actual reference bases at this position
            let bases = match read_ref_bases(fasta_path, &spec.chrom, position) {
                Ok(b) => b,
                Err(_) => {
                    // Fallback to random if reference read fails
                    (0..READ_LENGTH)
                        .map(|_| b"ACGT"[rng.random_range(0..4usize)])
                        .collect()
                }
            };

            reads.push(SynthRead {
                ref_id,
                position,
                bases,
                quality_scores: vec![30; READ_LENGTH as usize],
                mapq: 60,
                is_reverse: rng.random_range(0..2u8) == 1,
            });
        }
    }

    // Apply variant spikes
    for read in reads.iter_mut() {
        let read_start = read.position;
        let read_end = read_start + read.bases.len() as u32 - 1;

        let chrom_name = header
            .reference_sequences()
            .get_index(read.ref_id)
            .map(|(name, _)| std::str::from_utf8(name).unwrap_or(""))
            .unwrap_or("");

        for spike in spikes {
            if spike.chrom != chrom_name {
                continue;
            }
            if spike.pos < read_start || spike.pos > read_end {
                continue;
            }

            let offset = (spike.pos - read_start) as usize;
            if offset < read.bases.len() {
                let use_alt: bool = rng.random_range(0.0..1.0) < spike.allele_fraction;
                if use_alt {
                    read.bases[offset] = spike.alt_base;
                    // Apply spike overrides for alt-carrying reads
                    if let Some(strand) = spike.strand_override {
                        read.is_reverse = matches!(strand, Strand::Reverse);
                    }
                    if let Some(bq) = spike.alt_base_quality {
                        read.quality_scores[offset] = bq;
                    }
                    if let Some(mq) = spike.alt_mapq {
                        read.mapq = mq;
                    }
                }
                // If not using alt, base already has the reference base
            }
        }
    }

    reads.sort_by_key(|r| (r.ref_id, r.position));
    reads
}

/// Generate a BAI index for the BAM at `bam_path` using samtools.
/// Returns an error if samtools is not available.
pub fn index_bam(bam_path: &Path) -> io::Result<()> {
    let status = std::process::Command::new("samtools")
        .arg("index")
        .arg(bam_path)
        .status()?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "samtools index failed with {}",
            status
        )));
    }
    Ok(())
}

/// Fixture-sized region constants. Shared single source of truth for both
/// `fixture_specs()` here and `fixture_config()` in
/// `core/tests/fixture_helpers.rs`. If you change these, update both.
pub mod fixture_regions {
    /// PMP22 region  - same as production (35 kb).
    pub const PMP22: (u32, u32) = (15_229_777, 15_265_079);
    /// Chr2 control  - 100 kb (production: 20 MB).
    pub const CHR2_CONTROL: (u32, u32) = (50_000_000, 50_100_000);
    /// Boundary scan  - 550 kb (production: 3 MB).
    pub const BOUNDARY: (u32, u32) = (14_500_000, 15_050_000);
    /// CMT1A duplication start  - 400 kb dup region.
    pub const DUP_START: u32 = 14_600_000;
    /// CMT1A duplication end.
    pub const DUP_END: u32 = 15_000_000;
    /// Min classical run length for fixture-sized regions.
    pub const MIN_CLASSICAL_RUN_LENGTH: u32 = 300_000;
    /// Max classical run length for fixture-sized regions.
    pub const MAX_CLASSICAL_RUN_LENGTH: u32 = 600_000;
}

/// Returns the named fixture list for `cargo xtask synth-fixtures`.
/// Each entry is `(fixture_name, coverage_specs)`. Uses constants from
/// [`fixture_regions`] so `fixture_config()` can import the same values.
pub fn fixture_specs() -> Vec<(String, Vec<CoverageSpec>)> {
    use fixture_regions::*;

    let cs = |chrom: &str, start: u32, end: u32, depth: u32| CoverageSpec {
        chrom: chrom.to_string(),
        start,
        end,
        depth,
    };

    vec![
        (
            "cn2_30x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 30),
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 30),
                cs("chr17", BOUNDARY.0, BOUNDARY.1, 30),
            ],
        ),
        (
            "cn3_30x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 45),
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 30),
                cs("chr17", BOUNDARY.0, DUP_START - 1, 30),
                cs("chr17", DUP_START, DUP_END, 45),
                cs("chr17", DUP_END + 1, BOUNDARY.1, 30),
            ],
        ),
        (
            "cn1_30x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 15),
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 30),
                cs("chr17", BOUNDARY.0, DUP_START - 1, 30),
                cs("chr17", DUP_START, DUP_END, 15),
                cs("chr17", DUP_END + 1, BOUNDARY.1, 30),
            ],
        ),
        // --- Coverage-floor fixtures ---
        (
            "cn3_low_20x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 30), // PMP22 at dup depth
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 20), // 20x autosomal → Low
                cs("chr17", BOUNDARY.0, DUP_START - 1, 20),
                cs("chr17", DUP_START, DUP_END, 30),
                cs("chr17", DUP_END + 1, BOUNDARY.1, 20),
            ],
        ),
        (
            "cn3_floor_10x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 10),
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 10), // 10x autosomal → Refused
                cs("chr17", BOUNDARY.0, BOUNDARY.1, 10),
            ],
        ),
        (
            "cn2_26x".into(),
            vec![
                cs("chr17", PMP22.0, PMP22.1, 26),
                cs("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1, 26), // 26x autosomal → Full (25x boundary + margin for jitter)
                cs("chr17", BOUNDARY.0, BOUNDARY.1, 26),
            ],
        ),
    ]
}

/// High-level: generate a complete BAM + BAI from coverage specs.
pub fn generate_fixture(path: &Path, specs: &[CoverageSpec], seed: u64) -> io::Result<()> {
    let header = build_grch38_header();
    let reads = generate_reads(&header, specs, seed);
    write_bam(path, &header, &reads)?;
    index_bam(path)?;
    Ok(())
}

/// Generate all error-path fixtures into `out_dir`.
pub fn generate_error_fixtures(out_dir: &Path) -> io::Result<()> {
    use std::num::NonZeroUsize;

    std::fs::create_dir_all(out_dir)?;

    // 1. wrong_build.bam  - GRCh37 chr17 length in header.
    {
        let mut builder = noodles_sam::Header::builder();
        let len = NonZeroUsize::try_from(12_345_678usize).unwrap(); // wrong length (not GRCh37 or GRCh38)
        let ref_seq = Map::<ReferenceSequence>::new(len);
        builder = builder.add_reference_sequence("chr17", ref_seq);
        let header = builder.build();
        let specs = vec![CoverageSpec {
            chrom: "chr17".into(),
            start: 1_000_000,
            end: 1_010_000,
            depth: 5,
        }];
        let reads = generate_reads(&header, &specs, 42);
        write_bam(&out_dir.join("wrong_build.bam"), &header, &reads)?;
        // No BAI needed  - error fires before region query.
    }

    // 2. truncated.bam  - valid BAM with last 4 KB removed.
    {
        let header = build_grch38_header();
        let specs = vec![CoverageSpec {
            chrom: "chr17".into(),
            start: 15_229_777,
            end: 15_265_079,
            depth: 10,
        }];
        let reads = generate_reads(&header, &specs, 42);
        let path = out_dir.join("truncated.bam");
        write_bam(&path, &header, &reads)?;
        let len = std::fs::metadata(&path)?.len();
        let truncated_len = len.saturating_sub(4096);
        let file = std::fs::OpenOptions::new().write(true).open(&path)?;
        file.set_len(truncated_len)?;
    }

    // 3. empty_region.bam  - reads at chr1 only, none at PMP22/chr2.
    {
        let header = build_grch38_header();
        let specs = vec![CoverageSpec {
            chrom: "chr1".into(),
            start: 1_000_000,
            end: 1_010_000,
            depth: 10,
        }];
        let reads = generate_reads(&header, &specs, 42);
        write_bam(&out_dir.join("empty_region.bam"), &header, &reads)?;
        // No BAI needed.
    }

    // 4. missing_index/test.bam  - valid BAM in a subdirectory, no .bai.
    {
        let subdir = out_dir.join("missing_index");
        std::fs::create_dir_all(&subdir)?;
        let header = build_grch38_header();
        let specs = vec![CoverageSpec {
            chrom: "chr17".into(),
            start: 15_229_777,
            end: 15_265_079,
            depth: 5,
        }];
        let reads = generate_reads(&header, &specs, 42);
        write_bam(&subdir.join("test.bam"), &header, &reads)?;
        // Deliberately NO samtools index call.
    }

    // 5. invalid_format.bin  - not a BAM or CRAM.
    {
        let data: Vec<u8> = (0..1024).map(|i| (i % 256) as u8).collect();
        std::fs::write(out_dir.join("invalid_format.bin"), &data)?;
    }

    // 6. md5_mismatch.bam  - valid BAM with wrong M5 tag on chr17.
    {
        let mut builder = noodles_sam::Header::builder();
        let len = NonZeroUsize::try_from(nds_core::GRCH38_CHR17_LENGTH as usize).unwrap();
        let mut ref_seq = Map::<ReferenceSequence>::new(len);
        ref_seq.other_fields_mut().insert(
            noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM,
            "00000000000000000000000000000000".into(),
        );
        builder = builder.add_reference_sequence("chr17", ref_seq);
        let header = builder.build();
        let specs = vec![CoverageSpec {
            chrom: "chr17".into(),
            start: 15_229_777,
            end: 15_265_079,
            depth: 5,
        }];
        let reads = generate_reads(&header, &specs, 42);
        write_bam(&out_dir.join("md5_mismatch.bam"), &header, &reads)?;
    }

    Ok(())
}
