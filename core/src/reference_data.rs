//! Auditable reference data for every genomic coordinate used by
//! neuropathy-dna-scanner.
//!
//! Every coordinate, gene position, REP boundary, and disease association
//! carries the URL of its primary source in an `accessed` field. This module
//! is the single source of truth for the "magic numbers" in the tool.
//! hg19-to-hg38 conversions are cross-checked via UCSC liftOver.

use serde::Serialize;

/// A single genomic coordinate with its authoritative source.
#[derive(Debug, Clone, Serialize)]
pub struct SourcedCoordinate {
    pub value: u32,
    pub source: &'static str,
    pub source_url: &'static str,
    pub accessed: &'static str,
}

/// A genomic region with sourced start/end coordinates.
#[derive(Debug, Clone, Serialize)]
pub struct SourcedRegion {
    pub chrom: &'static str,
    pub start: SourcedCoordinate,
    pub end: SourcedCoordinate,
    pub description: &'static str,
}

/// A gene record with its source.
#[derive(Debug, Clone, Serialize)]
pub struct GeneRecord {
    pub symbol: &'static str,
    pub name: &'static str,
    pub region: SourcedRegion,
    pub omim: &'static str,
    pub omim_url: &'static str,
}

/// A REP (segmental duplication) element with its source.
#[derive(Debug, Clone, Serialize)]
pub struct RepRecord {
    pub name: &'static str,
    pub region: SourcedRegion,
    pub size_bp: u32,
    pub publication: &'static str,
    pub publication_url: &'static str,
}

/// A syndrome definition with OMIM reference.
#[derive(Debug, Clone, Serialize)]
pub struct SyndromeRecord {
    pub name: &'static str,
    pub short_name: &'static str,
    pub omim: &'static str,
    pub omim_url: &'static str,
    pub mechanism: &'static str,
    pub mediated_by: &'static str,
    pub expected_cn: u8,
    pub typical_size_description: &'static str,
    /// ICD-10 code. Source: WHO ICD-10 / Orphanet cross-reference.
    pub icd10: &'static str,
    /// ICD-11 code. Source: WHO ICD-11 / Orphanet cross-reference.
    pub icd11: &'static str,
    /// Orphanet disease number.
    pub orphanet: &'static str,
    pub orphanet_url: &'static str,
    /// GeneReviews URL (NCBI Bookshelf).
    pub genereviews_url: &'static str,
    /// GARD (Genetic and Rare Diseases Information Center) URL.
    pub gard_url: &'static str,
    /// One-paragraph description for non-medical readers.
    pub lay_description: &'static str,
    /// Common symptoms in plain language, sourced from GeneReviews or Orphanet.
    pub symptoms: &'static [&'static str],
    /// Citation for the symptom list.
    pub symptoms_source: &'static str,
    pub symptoms_source_url: &'static str,
}

/// Complete reference dataset for a given genome build.
#[derive(Debug, Clone, Serialize)]
pub struct ReferenceDataset {
    pub build: &'static str,
    pub build_source: &'static str,
    pub chr17_length: SourcedCoordinate,
    pub genes: Vec<GeneRecord>,
    pub reps: Vec<RepRecord>,
    pub syndromes: Vec<SyndromeRecord>,
    pub liftover_note: &'static str,
}

/// Syndrome records are build-independent (same diseases regardless of
/// reference genome). Shared between [`grch37`] and [`grch38`].
/// ICD codes verified via Orphanet cross-references; symptoms sourced from
/// GeneReviews (NCBI Bookshelf) and Orphanet.
fn syndromes() -> Vec<SyndromeRecord> {
    vec![
        SyndromeRecord {
            name: "Charcot-Marie-Tooth disease type 1A",
            short_name: "CMT1A",
            omim: "118220",
            omim_url: "https://omim.org/entry/118220",
            mechanism: "duplication",
            mediated_by: "CMT1A-REP",
            expected_cn: 3,
            typical_size_description: "~1.4 Mb duplication between CMT1A-REP elements",
            icd10: "G60.0",
            icd11: "LD41.G1",
            orphanet: "101081",
            orphanet_url: "https://www.orpha.net/en/disease/detail/101081",
            genereviews_url: "https://www.ncbi.nlm.nih.gov/books/NBK1358/",
            gard_url: "https://rarediseases.info.nih.gov/diseases/1245",
            lay_description: "CMT1A is the most common inherited nerve disease, affecting roughly 1 in 5,000 people. It is caused by having an extra copy of the PMP22 gene, which leads to progressive weakness and loss of sensation in the feet and hands. Most people notice symptoms in their teens or twenties, starting with difficulty walking, high-arched feet, and frequent ankle sprains. It progresses slowly over decades and does not affect life expectancy.",
            symptoms: &[
                "Progressive weakness in feet and lower legs",
                "High-arched feet (pes cavus)",
                "Difficulty walking, frequent tripping",
                "Loss of sensation in feet and hands",
                "Decreased reflexes",
                "Hammer toes",
                "Muscle wasting in lower legs",
            ],
            symptoms_source: "GeneReviews CMT Overview, revised Nov 2025",
            symptoms_source_url: "https://www.ncbi.nlm.nih.gov/books/NBK1358/",
        },
        SyndromeRecord {
            name: "Hereditary Neuropathy with liability to Pressure Palsies",
            short_name: "HNPP",
            omim: "162500",
            omim_url: "https://omim.org/entry/162500",
            mechanism: "deletion",
            mediated_by: "CMT1A-REP",
            expected_cn: 1,
            typical_size_description: "~1.4 Mb deletion (reciprocal of CMT1A)",
            icd10: "G60.0",
            icd11: "8C20.Y",
            orphanet: "640",
            orphanet_url: "https://www.orpha.net/en/disease/detail/640",
            genereviews_url: "https://www.ncbi.nlm.nih.gov/books/NBK1392/",
            gard_url: "https://rarediseases.info.nih.gov/diseases/6640",
            lay_description: "HNPP causes episodes of numbness and weakness when nerves are compressed, for example by leaning on an elbow or crossing legs. Most people recover fully within weeks, but repeated episodes can leave lasting mild weakness. It is caused by having only one copy of the PMP22 gene instead of the normal two.",
            symptoms: &[
                "Recurring episodes of numbness or tingling",
                "Painless muscle weakness triggered by pressure",
                "Foot drop (difficulty lifting the front of the foot)",
                "Carpal tunnel-like symptoms",
                "Episodes usually resolve within weeks",
            ],
            symptoms_source: "GeneReviews HNPP, revised 2023",
            symptoms_source_url: "https://www.ncbi.nlm.nih.gov/books/NBK1392/",
        },
        SyndromeRecord {
            name: "Potocki-Lupski syndrome",
            short_name: "PTLS",
            omim: "610883",
            omim_url: "https://omim.org/entry/610883",
            mechanism: "duplication",
            mediated_by: "SMS-REP",
            expected_cn: 3,
            typical_size_description: "~3.5 Mb duplication between SMS-REP elements",
            icd10: "Q92.3",
            icd11: "LD41.G1",
            orphanet: "1713",
            orphanet_url: "https://www.orpha.net/en/disease/detail/1713",
            genereviews_url: "https://www.ncbi.nlm.nih.gov/books/NBK447920/",
            gard_url: "https://rarediseases.info.nih.gov/diseases/10145",
            lay_description: "PTLS is a rare genetic condition caused by an extra copy of a region on chromosome 17 that includes the RAI1 gene. It typically causes developmental delays, especially in speech, mild to moderate intellectual disability, and sometimes heart defects. Many individuals have features of autism spectrum disorder.",
            symptoms: &[
                "Developmental delay, especially speech",
                "Mild to moderate intellectual disability",
                "Low muscle tone (hypotonia) in infancy",
                "Feeding difficulties in infancy",
                "Autism spectrum features (~60% of cases)",
                "Heart defects (~40% of cases)",
                "Sleep apnea",
            ],
            symptoms_source: "GeneReviews PTLS, revised 2023",
            symptoms_source_url: "https://www.ncbi.nlm.nih.gov/books/NBK447920/",
        },
        SyndromeRecord {
            name: "Smith-Magenis syndrome",
            short_name: "SMS",
            omim: "182290",
            omim_url: "https://omim.org/entry/182290",
            mechanism: "deletion",
            mediated_by: "SMS-REP",
            expected_cn: 1,
            typical_size_description: "~3.5 Mb deletion between SMS-REP elements",
            icd10: "Q93.5",
            icd11: "LD44.H1",
            orphanet: "819",
            orphanet_url: "https://www.orpha.net/en/disease/detail/819",
            genereviews_url: "https://www.ncbi.nlm.nih.gov/books/NBK1310/",
            gard_url: "https://rarediseases.info.nih.gov/diseases/7698",
            lay_description: "SMS is caused by a missing copy of the RAI1 gene region on chromosome 17. It causes intellectual disability, distinctive facial features, behavioral challenges (including self-injury and sleep disruption), and sometimes heart or kidney problems. The inverted sleep cycle (sleepy during the day, awake at night) is a hallmark feature.",
            symptoms: &[
                "Intellectual disability",
                "Distinctive facial features",
                "Inverted sleep-wake cycle",
                "Self-injurious behavior",
                "Speech and motor delays",
                "Short stature",
                "Heart defects (some cases)",
                "Kidney abnormalities (some cases)",
            ],
            symptoms_source: "Orphanet SMS overview",
            symptoms_source_url: "https://www.orpha.net/en/disease/detail/819",
        },
        SyndromeRecord {
            name: "Yuan-Harel-Lupski syndrome",
            short_name: "YUHAL",
            omim: "616652",
            omim_url: "https://omim.org/entry/616652",
            mechanism: "duplication",
            mediated_by: "CMT1A-REP + SMS-REP",
            expected_cn: 3,
            typical_size_description: "~6 Mb contiguous duplication spanning PMP22 through RAI1",
            icd10: "Q92.3",
            icd11: "",
            orphanet: "477817",
            orphanet_url: "https://www.orpha.net/en/disease/detail/477817",
            genereviews_url: "",
            gard_url: "https://rarediseases.info.nih.gov/diseases/17859",
            lay_description: "YUHAL is a very rare condition caused by a large duplication on chromosome 17 that includes both the PMP22 and RAI1 genes. It combines features of both CMT1A (peripheral nerve damage) and PTLS (developmental delays). The peripheral neuropathy may appear earlier and progress faster than in CMT1A alone.",
            symptoms: &[
                "Early-onset peripheral neuropathy (CMT1A features)",
                "Developmental delay (PTLS features)",
                "Intellectual disability",
                "Combined central and peripheral nervous system involvement",
            ],
            symptoms_source: "Orphanet YUHAL overview",
            symptoms_source_url: "https://www.orpha.net/en/disease/detail/477817",
        },
    ]
}

// ============================================================================
// Source constants
// ============================================================================

const NCBI_GENE_PMP22: &str = "https://www.ncbi.nlm.nih.gov/gene/5376";
const NCBI_GENE_RAI1: &str = "https://www.ncbi.nlm.nih.gov/gene/10743";
const NCBI_GRCH38_ASSEMBLY: &str = "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40";
const NCBI_GRCH37_ASSEMBLY: &str = "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25";
// UCSC Segmental Duplications track (Eichler lab WGAC pipeline).
// This is the actual database that contains the CMT1A-REP and SMS-REP
// coordinates. Inoue 2001 / Park 2002 describe the biology but the
// papers themselves do not contain base-pair coordinates.
const UCSC_SEGDUPS_HG38: &str = "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=range&position=chr17:14100000-15700000&hgta_outputType=primaryTable&hgta_doTopSubmit=submit";
const UCSC_SEGDUPS_HG19: &str = "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=range&position=chr17:14000000-15600000&hgta_outputType=primaryTable&hgta_doTopSubmit=submit";
const UCSC_SEGDUPS_SMS_HG38: &str = "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=range&position=chr17:16600000-20600000&hgta_outputType=primaryTable&hgta_doTopSubmit=submit";
const UCSC_SEGDUPS_SMS_HG19: &str = "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_regionType=range&position=chr17:16500000-20400000&hgta_outputType=primaryTable&hgta_doTopSubmit=submit";
const ACCESSED_DATE: &str = "2026-04-20";

// ============================================================================
// GRCh38 (hg38) dataset
// ============================================================================

/// GRCh38 reference dataset with all coordinates verified from primary sources.
///
/// Gene positions: NCBI Gene database, GRCh38.p14 annotation.
/// CMT1A-REP positions: UCSC genomicSuperDups track (hg38), Eichler lab WGAC.
/// SMS-REP positions: UCSC genomicSuperDups track (hg38) for >100 kb entries.
/// SMS region boundaries cross-checked with ClinVar RCV000591005 (liftover).
pub fn grch38() -> ReferenceDataset {
    ReferenceDataset {
        build: "GRCh38",
        build_source: NCBI_GRCH38_ASSEMBLY,
        chr17_length: SourcedCoordinate {
            value: 83_257_441,
            source: "NCBI Assembly GRCh38.p14 sequence report",
            source_url: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
            accessed: ACCESSED_DATE,
        },
        genes: vec![
            GeneRecord {
                symbol: "PMP22",
                name: "Peripheral Myelin Protein 22",
                region: SourcedRegion {
                    chrom: "chr17",
                    start: SourcedCoordinate { value: 15_229_779, source: "NCBI Gene 5376, GRCh38.p14", source_url: NCBI_GENE_PMP22, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 15_265_326, source: "NCBI Gene 5376, GRCh38.p14", source_url: NCBI_GENE_PMP22, accessed: ACCESSED_DATE },
                    description: "PMP22 gene locus (complement strand)",
                },
                omim: "601097",
                omim_url: "https://omim.org/entry/601097",
            },
            GeneRecord {
                symbol: "RAI1",
                name: "Retinoic Acid Induced 1",
                region: SourcedRegion {
                    chrom: "chr17",
                    start: SourcedCoordinate { value: 17_681_458, source: "NCBI Gene 10743, GRCh38.p14", source_url: NCBI_GENE_RAI1, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 17_811_453, source: "NCBI Gene 10743, GRCh38.p14", source_url: NCBI_GENE_RAI1, accessed: ACCESSED_DATE },
                    description: "RAI1 gene locus",
                },
                omim: "607642",
                omim_url: "https://omim.org/entry/607642",
            },
        ],
        reps: vec![
            RepRecord {
                name: "CMT1A-REP distal (telomeric)",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg38: 0-based 14170711-14194597 = 1-based 14170712-14194597
                    start: SourcedCoordinate { value: 14_170_712, source: "UCSC genomicSuperDups hg38 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG38, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 14_194_597, source: "UCSC genomicSuperDups hg38 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG38, accessed: ACCESSED_DATE },
                    description: "Distal (telomeric) CMT1A-REP element - NAHR substrate for CMT1A/HNPP",
                },
                size_bp: 23_886,
                publication: "UCSC genomicSuperDups hg38",
                publication_url: UCSC_SEGDUPS_HG38,
            },
            RepRecord {
                name: "CMT1A-REP proximal (centromeric)",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg38: 0-based 15567588-15591460 = 1-based 15567589-15591460
                    start: SourcedCoordinate { value: 15_567_589, source: "UCSC genomicSuperDups hg38 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG38, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 15_591_460, source: "UCSC genomicSuperDups hg38 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG38, accessed: ACCESSED_DATE },
                    description: "Proximal (centromeric) CMT1A-REP element - NAHR substrate for CMT1A/HNPP",
                },
                size_bp: 23_872,
                publication: "UCSC genomicSuperDups hg38",
                publication_url: UCSC_SEGDUPS_HG38,
            },
            RepRecord {
                name: "SMS-REP distal",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg38: 0-based 16685735-16813672 = 1-based 16685736-16813672
                    start: SourcedCoordinate { value: 16_685_736, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 16_813_672, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    description: "Distal SMS-REP element (~128 kb) - NAHR substrate for PTLS/SMS",
                },
                size_bp: 127_937,
                publication: "UCSC genomicSuperDups hg38",
                publication_url: UCSC_SEGDUPS_SMS_HG38,
            },
            RepRecord {
                name: "SMS-REP middle",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg38: 0-based 18384909-18568180 = 1-based 18384910-18568180
                    start: SourcedCoordinate { value: 18_384_910, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 18_568_180, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    description: "Middle SMS-REP element (~183 kb) - involved in atypical SMS deletions",
                },
                size_bp: 183_271,
                publication: "UCSC genomicSuperDups hg38",
                publication_url: UCSC_SEGDUPS_SMS_HG38,
            },
            RepRecord {
                name: "SMS-REP proximal",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg38: paired target region 20316336-20560909
                    // Using the target (proximal) coordinates from the distal-proximal alignment
                    start: SourcedCoordinate { value: 20_316_337, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, paired target of distal REP)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 20_560_909, source: "UCSC genomicSuperDups hg38 (Eichler WGAC, paired target of middle REP)", source_url: UCSC_SEGDUPS_SMS_HG38, accessed: ACCESSED_DATE },
                    description: "Proximal SMS-REP element (~245 kb) - NAHR substrate for PTLS/SMS",
                },
                size_bp: 244_573,
                publication: "UCSC genomicSuperDups hg38",
                publication_url: UCSC_SEGDUPS_SMS_HG38,
            },
        ],
        syndromes: syndromes(),
        liftover_note: "Gene positions from NCBI Gene (GRCh38.p14). REP coordinates from UCSC Segmental Duplications track (genomicSuperDups, hg38 assembly, Eichler lab WGAC pipeline). Each coordinate links directly to the UCSC Table Browser query that produces it.",
    }
}

// ============================================================================
// GRCh37 (hg19) dataset
// ============================================================================

/// GRCh37 reference dataset with verified coordinates.
///
/// Gene positions: NCBI Gene database, GRCh37.p13 annotation.
/// CMT1A-REP positions: UCSC genomicSuperDups track (hg19), Eichler lab WGAC.
/// SMS-REP positions: UCSC genomicSuperDups track (hg19) for >100 kb entries;
/// overall PTLS region boundary from ClinVar RCV000591005.
pub fn grch37() -> ReferenceDataset {
    ReferenceDataset {
        build: "GRCh37",
        build_source: NCBI_GRCH37_ASSEMBLY,
        chr17_length: SourcedCoordinate {
            value: 81_195_210,
            source: "NCBI Assembly GRCh37.p13 sequence report",
            source_url: "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt",
            accessed: ACCESSED_DATE,
        },
        genes: vec![
            GeneRecord {
                symbol: "PMP22",
                name: "Peripheral Myelin Protein 22",
                region: SourcedRegion {
                    chrom: "chr17",
                    start: SourcedCoordinate { value: 15_133_096, source: "NCBI Gene 5376, GRCh37.p13", source_url: NCBI_GENE_PMP22, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 15_168_643, source: "NCBI Gene 5376, GRCh37.p13", source_url: NCBI_GENE_PMP22, accessed: ACCESSED_DATE },
                    description: "PMP22 gene locus (complement strand)",
                },
                omim: "601097",
                omim_url: "https://omim.org/entry/601097",
            },
            GeneRecord {
                symbol: "RAI1",
                name: "Retinoic Acid Induced 1",
                region: SourcedRegion {
                    chrom: "chr17",
                    start: SourcedCoordinate { value: 17_584_772, source: "NCBI Gene 10743, GRCh37.p13", source_url: NCBI_GENE_RAI1, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 17_714_767, source: "NCBI Gene 10743, GRCh37.p13", source_url: NCBI_GENE_RAI1, accessed: ACCESSED_DATE },
                    description: "RAI1 gene locus",
                },
                omim: "607642",
                omim_url: "https://omim.org/entry/607642",
            },
        ],
        reps: vec![
            RepRecord {
                name: "CMT1A-REP distal (telomeric)",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg19: 0-based 14074028-14097914 = 1-based 14074029-14097914
                    start: SourcedCoordinate { value: 14_074_029, source: "UCSC genomicSuperDups hg19 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG19, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 14_097_914, source: "UCSC genomicSuperDups hg19 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG19, accessed: ACCESSED_DATE },
                    description: "Distal (telomeric) CMT1A-REP element",
                },
                size_bp: 23_886,
                publication: "UCSC genomicSuperDups hg19",
                publication_url: UCSC_SEGDUPS_HG19,
            },
            RepRecord {
                name: "CMT1A-REP proximal (centromeric)",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg19: 0-based 15470902-15494774 = 1-based 15470903-15494774
                    start: SourcedCoordinate { value: 15_470_903, source: "UCSC genomicSuperDups hg19 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG19, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 15_494_774, source: "UCSC genomicSuperDups hg19 (Eichler WGAC)", source_url: UCSC_SEGDUPS_HG19, accessed: ACCESSED_DATE },
                    description: "Proximal (centromeric) CMT1A-REP element",
                },
                size_bp: 23_872,
                publication: "UCSC genomicSuperDups hg19",
                publication_url: UCSC_SEGDUPS_HG19,
            },
            RepRecord {
                name: "SMS-REP distal",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg19: 0-based 16589049-16716986 = 1-based 16589050-16716986
                    start: SourcedCoordinate { value: 16_589_050, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 16_716_986, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    description: "Distal SMS-REP element (~128 kb) - NAHR substrate for PTLS/SMS",
                },
                size_bp: 127_937,
                publication: "UCSC genomicSuperDups hg19",
                publication_url: UCSC_SEGDUPS_SMS_HG19,
            },
            RepRecord {
                name: "SMS-REP middle",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg19: 0-based 18288223-18471494 = 1-based 18288224-18471494
                    start: SourcedCoordinate { value: 18_288_224, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 18_471_494, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, >100 kb entry)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    description: "Middle SMS-REP element (~183 kb) - involved in atypical SMS deletions",
                },
                size_bp: 183_271,
                publication: "UCSC genomicSuperDups hg19",
                publication_url: UCSC_SEGDUPS_SMS_HG19,
            },
            RepRecord {
                name: "SMS-REP proximal",
                region: SourcedRegion {
                    chrom: "chr17",
                    // UCSC genomicSuperDups hg19: paired targets from distal+middle REP alignments
                    // 0-based 20219649-20464222 = 1-based 20219650-20464222
                    start: SourcedCoordinate { value: 20_219_650, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, paired target of distal REP)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    end: SourcedCoordinate { value: 20_464_222, source: "UCSC genomicSuperDups hg19 (Eichler WGAC, paired target of middle REP)", source_url: UCSC_SEGDUPS_SMS_HG19, accessed: ACCESSED_DATE },
                    description: "Proximal SMS-REP element (~245 kb) - NAHR substrate for PTLS/SMS",
                },
                size_bp: 244_573,
                publication: "UCSC genomicSuperDups hg19",
                publication_url: UCSC_SEGDUPS_SMS_HG19,
            },
        ],
        syndromes: syndromes(),
        liftover_note: "Gene positions from NCBI Gene (GRCh37.p13). REP coordinates from UCSC Segmental Duplications track (genomicSuperDups, hg19 assembly, Eichler lab WGAC pipeline). Each coordinate links directly to the UCSC Table Browser query that produces it.",
    }
}

/// Select the appropriate reference dataset based on the detected chr17 length.
pub fn for_chr17_length(length: u64) -> Option<ReferenceDataset> {
    if length == 83_257_441 {
        Some(grch38())
    } else if length == 81_195_210 {
        Some(grch37())
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grch38_source_urls_are_valid() {
        let ds = grch38();
        assert!(ds.chr17_length.source_url.starts_with("https://"));
        for gene in &ds.genes {
            assert!(
                gene.region.start.source_url.starts_with("https://"),
                "gene {} start",
                gene.symbol
            );
            assert!(
                gene.region.end.source_url.starts_with("https://"),
                "gene {} end",
                gene.symbol
            );
            assert!(gene.omim_url.starts_with("https://"));
        }
        for rep in &ds.reps {
            assert!(
                rep.region.start.source_url.starts_with("https://"),
                "rep {} start",
                rep.name
            );
            assert!(
                rep.region.end.source_url.starts_with("https://"),
                "rep {} end",
                rep.name
            );
            assert!(rep.publication_url.starts_with("https://"));
        }
        for syn in &ds.syndromes {
            assert!(
                syn.omim_url.starts_with("https://"),
                "syndrome {} omim",
                syn.name
            );
        }
    }

    #[test]
    fn grch37_source_urls_are_valid() {
        let ds = grch37();
        assert!(ds.chr17_length.source_url.starts_with("https://"));
        for gene in &ds.genes {
            assert!(gene.region.start.source_url.starts_with("https://"));
        }
        for rep in &ds.reps {
            assert!(rep.region.start.source_url.starts_with("https://"));
        }
    }

    #[test]
    fn chr17_lengths_match_known_values() {
        assert_eq!(grch38().chr17_length.value, 83_257_441);
        assert_eq!(grch37().chr17_length.value, 81_195_210);
    }

    #[test]
    fn for_chr17_length_selects_correct_build() {
        assert_eq!(for_chr17_length(83_257_441).unwrap().build, "GRCh38");
        assert_eq!(for_chr17_length(81_195_210).unwrap().build, "GRCh37");
        assert!(for_chr17_length(12345).is_none());
    }

    #[test]
    fn pmp22_is_inside_cmt1a_reps() {
        let ds = grch38();
        let pmp22 = &ds.genes[0];
        let distal_rep = &ds.reps[0];
        let proximal_rep = &ds.reps[1];
        assert_eq!(pmp22.symbol, "PMP22");
        // PMP22 gene should be between the two CMT1A-REP elements
        assert!(
            pmp22.region.start.value > distal_rep.region.end.value,
            "PMP22 start should be past distal REP end"
        );
        assert!(
            pmp22.region.end.value < proximal_rep.region.start.value,
            "PMP22 end should be before proximal REP start"
        );
    }

    #[test]
    fn rai1_is_inside_sms_reps() {
        let ds = grch38();
        let rai1 = &ds.genes[1];
        let sms_distal = &ds.reps[2];
        let sms_proximal = &ds.reps[4];
        assert_eq!(rai1.symbol, "RAI1");
        // RAI1 should be between SMS-REP distal and proximal
        assert!(
            rai1.region.start.value > sms_distal.region.start.value,
            "RAI1 start should be past SMS-REP distal start"
        );
        assert!(
            rai1.region.end.value < sms_proximal.region.end.value,
            "RAI1 end should be before SMS-REP proximal end"
        );
    }

    #[test]
    fn five_syndromes_defined() {
        assert_eq!(grch38().syndromes.len(), 5);
        assert_eq!(grch37().syndromes.len(), 5);
    }

    #[test]
    fn datasets_have_same_structure() {
        let g38 = grch38();
        let g37 = grch37();
        assert_eq!(g38.genes.len(), g37.genes.len());
        assert_eq!(g38.reps.len(), g37.reps.len());
        assert_eq!(g38.syndromes.len(), g37.syndromes.len());
        // Same gene symbols
        for (a, b) in g38.genes.iter().zip(g37.genes.iter()) {
            assert_eq!(a.symbol, b.symbol);
        }
    }
}
