/** TypeScript types for the WASM boundary. */

export interface WasmModule {
  version(): string;
  render_matrix_text(): string;
  render_matrix_html(): string;
  detect_sex_from_index(crai_bytes: Uint8Array, header_bytes: Uint8Array): string;
  required_ranges_bam(
    bai_bytes: Uint8Array,
    header_bytes: Uint8Array,
    regions_json: string,
  ): RangesReply;
  required_ranges_cram(
    crai_bytes: Uint8Array,
    header_bytes: Uint8Array,
    regions_json: string,
  ): RangesReply;
  analyze_bam(
    bam_bytes: Uint8Array,
    bai_bytes: Uint8Array,
    config_json: string,
  ): AnalyzeReply;
  analyze_bam_sparse(
    data: Uint8Array,
    ranges_json: string,
    index_bytes: Uint8Array,
    config_json: string,
  ): AnalyzeReply;
  analyze_cram_sparse(
    data: Uint8Array,
    ranges_json: string,
    index_bytes: Uint8Array,
    ref_chr17: Uint8Array,
    ref_chr2: Uint8Array,
    config_json: string,
  ): AnalyzeReply;
  detect_build(header_bytes: Uint8Array, is_cram: boolean): string;
  reference_data_json(build: string): string;
  screening_manifest(build: string): string;
  catalog_csv(): string;
  catalog_info(): string;
  catalog_regions(): string;
  required_ranges_cram_tight(
    crai_bytes: Uint8Array,
    header_bytes: Uint8Array,
    regions_json: string,
  ): RangesReply;
  check_variants_bam(bam_bytes: Uint8Array, bai_bytes: Uint8Array): PileupResult;
  check_variants_bam_sparse(
    data: Uint8Array,
    ranges_json: string,
    bai_bytes: Uint8Array,
  ): PileupResult;
  check_variants_cram(
    cram_bytes: Uint8Array,
    crai_bytes: Uint8Array,
    ref_chr17: Uint8Array,
    ref_chr2: Uint8Array,
    sex: string,
  ): PileupResult;
  check_variants_cram_chrom(
    cram_data: Uint8Array,
    ranges_json: string,
    crai_bytes: Uint8Array,
    ref_seq: Uint8Array,
    chrom_name: string,
    sex: string,
  ): PileupResult;
  check_variants_cram_sparse(
    data: Uint8Array,
    ranges_json: string,
    crai_bytes: Uint8Array,
    ref_chr17: Uint8Array,
    ref_chr2: Uint8Array,
    sex: string,
  ): PileupResult;
}

export interface RangesReply {
  status: "ok" | "error";
  ranges?: ByteRange[];
  error?: string;
}

export interface ByteRange {
  offset: number;
  length: number;
}

export interface AnalyzeReply {
  status: "ok" | "error";
  depth?: DepthResult;
  boundary?: BoundaryResult;
  interpretation?: InterpretationResult;
  progress?: ProgressEvent[];
  header_sequences?: number;
  elapsed_ms?: number;
  message?: string;
}

export interface DepthResult {
  pmp22_mean_depth: number;
  autosomal_mean_depth: number;
  ratio: number;
  estimated_cn: number;
  rai1_mean_depth?: number;
  rai1_ratio?: number;
  rai1_estimated_cn?: number;
}

export interface BoundaryResult {
  window_count: number;
  duplicated_windows: number;
  deleted_windows: number;
  runs: BoundaryRun[];
}

export interface BoundaryRun {
  start: number;
  end: number;
  length: number;
  direction: string;
  mean_normalized_depth: number;
  bridged_gaps: number;
}

export interface InterpretationResult {
  copy_number: number | null;
  confidence: string;
  subtype_call: string | null;
  plain_language: string;
}

export interface ProgressEvent {
  kind: string;
  phase?: string;
  region?: string;
  records_seen?: number;
  last_position?: number;
  pmp22_mean_depth?: number;
  autosomal_mean_depth?: number;
  ratio?: number;
  estimated_cn?: number;
  windows?: number;
  longest_run_length?: number | null;
  confidence?: string;
  subtype_call?: string | null;
}

export interface RegionInput {
  chrom: string;
  start: number;
  end: number;
}

/** Reference data JSON shape from reference_data_json(). */
export interface ReferenceDataset {
  build: string;
  build_source: string;
  genes: GeneRecord[];
  reps: RepRecord[];
  syndromes: SyndromeRecord[];
  liftover_note?: string;
}

export interface SourcedCoordinate {
  value: number;
  source: string;
  source_url: string;
}

export interface GeneRecord {
  symbol: string;
  name: string;
  region: {
    chrom: string;
    start: SourcedCoordinate;
    end: SourcedCoordinate;
  };
}

export interface RepRecord {
  name: string;
  region: {
    chrom: string;
    start: SourcedCoordinate;
    end: SourcedCoordinate;
  };
  size_bp: number;
  publication: string;
  publication_url: string;
}

/** Pileup genotyping result from check_variants_bam(). */
export interface PileupResult {
  positive_findings: CheckedVariant[];
  total_checked: number;
  no_coverage_count: number;
  low_coverage_count: number;
  catalog_date: string;
  catalog_gene_count: number;
  catalog_variant_count: number;
}

export interface CheckedVariant {
  gene: string;
  chrom: string;
  pos: number;
  ref_allele: string;
  alt_allele: string;
  clinvar_id: number;
  significance: string;
  review_stars: number;
  condition: string;
  inheritance: string;
  consequence: string;
  call: VariantCall;
}

/**
 * Mirror of `nds_core::pileup::VariantCall` after Serde serialisation.
 * Serde serialises enum variants with data as `{ VariantName: payload }`,
 * and unit variants as `"VariantName"`. Keep in sync with pileup.rs.
 */
interface Counts {
  ref_count: number;
  alt_count: number;
  allele_balance: number;
}

export type VariantCall =
  | "NoCoverage"
  | { ReferenceOnly: { depth: number } }
  | { Heterozygous: Counts }
  | { Homozygous: Counts }
  | { Marginal: Counts }
  | {
      StrandBias: Counts & { forward_alt: number; reverse_alt: number };
    }
  | { LowCoverage: { depth: number } }
  | { Hemizygous: { ref_count: number; alt_count: number } };

export interface ScreeningManifest {
  sv_screening: SvScreeningItem[];
  cmt_core_genes: SnvScreeningItem[];
  extended_panel_genes: SnvScreeningItem[];
  catalog_date: string;
  catalog_variant_count: number;
  catalog_snv_count: number;
  catalog_skipped_indels: number;
  catalog_gene_count: number;
  catalog_filter: string;
  catalog_source: string;
}

export interface SvScreeningItem {
  type: "sv";
  short_name: string;
  name: string;
  gene: string;
  mechanism: string;
  omim: string;
  omim_url: string;
  icd10?: string;
  icd11?: string;
  orphanet?: string;
  orphanet_url?: string;
  genereviews_url?: string;
  extra_check?: boolean;
}

export interface SnvScreeningItem {
  gene: string;
  chrom: string;
  condition: string;
  variant_count: number;
  inheritance: string;
  is_cmt_core: boolean;
}

export interface CatalogMeta {
  generated_date: string;
  clinvar_date: string;
  panelapp_panel: string;
  gene_count: number;
  variant_count: number;
  min_review_stars: number;
}

export interface SyndromeRecord {
  short_name: string;
  name: string;
  gene: string;
  mechanism: string;
  expected_cn: number;
  mediated_by: string;
  omim: string;
  omim_url: string;
  icd10?: string;
  icd11?: string;
  orphanet?: string;
  orphanet_url?: string;
  genereviews_url?: string;
  gard_url?: string;
  lay_description?: string;
  symptoms?: string[];
  symptoms_source?: string;
  symptoms_source_url?: string;
}
