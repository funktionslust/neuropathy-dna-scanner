import { ref, computed, reactive } from "vue";
import type {
  WasmModule,
  AnalyzeReply,
  CheckedVariant,
  PileupResult,
  RegionInput,
  ProgressEvent,
} from "../types/wasm";

export type AnalysisState = "idle" | "loading" | "analyzing" | "done" | "error";

export interface DetectedFiles {
  bam: File | null;
  bai: File | null;
  cram: File | null;
  crai: File | null;
  ref: File | null;
}

export interface ProgressInfo {
  phase: string;
  detail: string;
  percent: number; // 0-100, -1 = indeterminate
  step: number; // current step (1-based)
  totalSteps: number; // total steps
}

export type ChromStatus = "pending" | "loading" | "done" | "skipped" | "error";

export interface ChromScreeningState {
  chrom: string;
  status: ChromStatus;
  variantCount: number;
  findings: number;
}

export interface ScreeningProgress {
  chromosomes: ChromScreeningState[];
  currentChrom: string;
}

const state = ref<AnalysisState>("idle");
const logs = ref<string[]>([]);
const result = ref<AnalyzeReply | null>(null);
const pileupResult = ref<PileupResult | null>(null);
const screeningProgress = reactive<ScreeningProgress>({
  chromosomes: [],
  currentChrom: "",
});
const detectedBuild = ref<string | null>(null);
const detectedSexValue = ref<string>("unknown");
const progress = reactive<ProgressInfo>({
  phase: "",
  detail: "",
  percent: 0,
  step: 0,
  totalSteps: 0,
});

// Cached reference sequences for CRAM
let cachedRefChr17: Uint8Array | null = null;
let cachedRefChr2: Uint8Array | null = null;

function log(msg: string) {
  logs.value.push(msg);
}

function setProgress(
  phase: string,
  detail: string,
  percent: number,
  step?: number,
  totalSteps?: number,
) {
  progress.phase = phase;
  progress.detail = detail;
  progress.percent = percent;
  if (step !== undefined) progress.step = step;
  if (totalSteps !== undefined) progress.totalSteps = totalSteps;
}

export function formatSize(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1048576) return `${(bytes / 1024).toFixed(1)} KB`;
  if (bytes < 1073741824) return `${(bytes / 1048576).toFixed(1)} MB`;
  return `${(bytes / 1073741824).toFixed(2)} GB`;
}

// Yield to the browser event loop so the UI can repaint.
// Double-yield via requestAnimationFrame + setTimeout ensures the
// browser actually renders before we resume.
function yieldToUI(): Promise<void> {
  return new Promise((resolve) => {
    requestAnimationFrame(() => setTimeout(resolve, 0));
  });
}

// ---------------------------------------------------------------------------
// FASTA reference loading - one-pass index + windowed loading
// ---------------------------------------------------------------------------

/** Cached FASTA index: chrom -> byte offsets in the file. */
let fastaIndex: Map<string, { seqStart: number; seqEnd: number }> | null =
  null;

/**
 * Build an index of ALL chromosome offsets in ONE pass through the FASTA.
 * ~30s for 3 GB but only done once. After this, any chromosome or
 * sub-region can be loaded instantly by byte offset.
 */
async function buildFastaIndex(
  file: File,
): Promise<Map<string, { seqStart: number; seqEnd: number }>> {
  if (fastaIndex) return fastaIndex;

  log("Building FASTA index (one-time scan)...");
  const chunkSize = 4 * 1024 * 1024;
  const headers: { name: string; seqStart: number }[] = [];

  // Check first line
  const first = await file.slice(0, 4096).arrayBuffer();
  const firstText = new TextDecoder("ascii").decode(new Uint8Array(first));
  if (firstText.startsWith(">")) {
    const nl = firstText.indexOf("\n");
    const sp = firstText.indexOf(" ");
    const nameEnd = Math.min(
      sp >= 0 ? sp : firstText.length,
      nl >= 0 ? nl : firstText.length,
    );
    const name = firstText.substring(1, nameEnd);
    if (!name.includes("_")) {
      headers.push({ name, seqStart: nl + 1 });
    }
  }

  // Scan for all ">" headers
  for (let offset = 0; offset < file.size; offset += chunkSize - 200) {
    if (offset % (64 * 1024 * 1024) < chunkSize) {
      setProgress(
        "Indexing reference genome",
        `${formatSize(offset)} / ${formatSize(file.size)}`,
        -1,
      );
    }
    const end = Math.min(offset + chunkSize, file.size);
    const buf = await file.slice(offset, end).arrayBuffer();
    const text = new TextDecoder("ascii", { fatal: false }).decode(
      new Uint8Array(buf),
    );

    let searchFrom = 0;
    while (true) {
      const idx = text.indexOf("\n>", searchFrom);
      if (idx < 0) break;
      const hStart = idx + 2;
      const sp = text.indexOf(" ", hStart);
      const nl = text.indexOf("\n", hStart);
      const nameEnd = Math.min(
        sp >= 0 ? sp : text.length,
        nl >= 0 ? nl : text.length,
      );
      const name = text.substring(hStart, nameEnd);
      const seqStart = offset + (nl >= 0 ? nl + 1 : nameEnd);
      if (!name.includes("_") && !headers.some((h) => h.name === name)) {
        headers.push({ name, seqStart });
      }
      searchFrom = hStart + 1;
    }
  }

  // Build index: each chrom ends where the next begins
  const index = new Map<string, { seqStart: number; seqEnd: number }>();
  headers.sort((a, b) => a.seqStart - b.seqStart);
  for (let i = 0; i < headers.length; i++) {
    const h = headers[i]!;
    // seqEnd = byte before the next header line (the "\n>" marker)
    const seqEnd =
      i + 1 < headers.length ? headers[i + 1]!.seqStart : file.size;
    index.set(h.name, { seqStart: h.seqStart, seqEnd });
  }

  fastaIndex = index;
  log(`FASTA index: ${index.size} chromosomes`);
  return index;
}

/**
 * Look up a chromosome's byte range in the FASTA.
 * Uses the cached index (builds it on first call).
 */
async function findChromInFasta(
  file: File,
  chromName: string,
): Promise<{ seqStart: number; seqEnd: number }> {
  const index = await buildFastaIndex(file);
  const entry = index.get(chromName);
  if (!entry) throw new Error(`${chromName} not found in FASTA`);
  return entry;
}

/**
 * Load only a WINDOW of reference sequence from the FASTA.
 * `genomeStart` and `genomeEnd` are 0-based genomic coordinates.
 * The FASTA has ~60 bases per line + newline, so we compute byte
 * offsets from the chromosome's seqStart.
 */
async function readRefWindow(
  file: File,
  chromName: string,
  genomeStart: number,
  genomeEnd: number,
): Promise<Uint8Array> {
  const index = await buildFastaIndex(file);
  const entry = index.get(chromName);
  if (!entry) throw new Error(`${chromName} not in FASTA index`);

  // FASTA format: 60 bases per line, 61 bytes per line (60 + newline).
  // Byte offset for genomic position P within the chromosome:
  //   byteOffset = seqStart + P + floor(P / 60)
  // (each 60 bases has one extra newline byte)
  const bytesPerLine = 61; // 60 bases + \n
  const basesPerLine = 60;

  const byteStart =
    entry.seqStart +
    genomeStart +
    Math.floor(genomeStart / basesPerLine);
  const byteEnd = Math.min(
    entry.seqStart +
      genomeEnd +
      Math.floor(genomeEnd / basesPerLine) +
      2,
    entry.seqEnd,
  );

  const buf = await file.slice(byteStart, byteEnd).arrayBuffer();
  const bytes = new Uint8Array(buf);

  // Strip newlines
  const result = new Uint8Array(genomeEnd - genomeStart + 1);
  let writePos = 0;
  for (let i = 0; i < bytes.length && writePos < result.length; i++) {
    const b = bytes[i]!;
    if (b !== 10 && b !== 13) {
      result[writePos++] = b;
    }
  }
  return result.subarray(0, writePos);
}

async function readChromSequence(
  file: File,
  chromName: string,
  seqStart: number,
  seqEnd: number,
): Promise<Uint8Array> {
  const rawSize = seqEnd - seqStart;

  const estimatedSize = Math.ceil(rawSize * 0.985);
  const resultBuf = new Uint8Array(estimatedSize);
  let writePos = 0;
  // Use 8 MB chunks for more frequent UI yields (32 MB blocks the
  // browser for ~2-3s, triggering "page unresponsive" dialogs).
  const chunkSize = 8 * 1024 * 1024;

  for (let offset = seqStart; offset < seqEnd; offset += chunkSize) {
    const end = Math.min(offset + chunkSize, seqEnd);
    const buf = await file.slice(offset, end).arrayBuffer();
    const bytes = new Uint8Array(buf);
    for (let i = 0; i < bytes.length; i++) {
      const b = bytes[i]!;
      if (b !== 10 && b !== 13) {
        resultBuf[writePos++] = b;
      }
    }
    // Yield to UI but don't override the caller's progress label
    await yieldToUI();
  }

  return resultBuf.subarray(0, writePos);
}

async function loadReferenceFromFasta(fastaFile: File): Promise<void> {
  if (cachedRefChr17 && cachedRefChr2) {
    log("Reference loaded from cache.");
    return;
  }

  // CRAM always has 6 steps; reference loading is step 4
  setProgress("Loading reference genome", `Scanning ${fastaFile.name}...`, 55, 4, 6);
  log(`Loading reference from ${fastaFile.name} (${formatSize(fastaFile.size)})`);

  const chr17Pos = await findChromInFasta(fastaFile, "chr17");
  const chr2Pos = await findChromInFasta(fastaFile, "chr2");

  cachedRefChr17 = await readChromSequence(
    fastaFile,
    "chr17",
    chr17Pos.seqStart,
    chr17Pos.seqEnd,
  );
  cachedRefChr2 = await readChromSequence(
    fastaFile,
    "chr2",
    chr2Pos.seqStart,
    chr2Pos.seqEnd,
  );

  log(
    `Reference loaded: chr17=${formatSize(cachedRefChr17.byteLength)}, chr2=${formatSize(cachedRefChr2.byteLength)}`,
  );
}

// ---------------------------------------------------------------------------
// File.slice() with retry
// ---------------------------------------------------------------------------

async function sliceWithRetry(
  file: File,
  offset: number,
  length: number,
  retries = 3,
): Promise<ArrayBuffer> {
  for (let attempt = 1; attempt <= retries; attempt++) {
    try {
      await new Promise((resolve) => setTimeout(resolve, 50));
      return await file.slice(offset, offset + length).arrayBuffer();
    } catch (e) {
      if (attempt < retries) {
        log(
          `Retry ${attempt}/${retries} for offset ${offset.toLocaleString()}`,
        );
        await new Promise((resolve) => setTimeout(resolve, 500 * attempt));
      } else {
        throw e;
      }
    }
  }
  throw new Error("unreachable");
}

// ---------------------------------------------------------------------------
// Main analysis pipeline
// ---------------------------------------------------------------------------

async function runAnalysis(
  wasm: WasmModule,
  detected: DetectedFiles,
): Promise<void> {
  const isCram = !!(detected.cram && detected.crai);
  const isBam = !!(detected.bam && detected.bai);
  if (!isCram && !isBam) throw new Error("BAM+BAI or CRAM+CRAI required");

  const dataFile = isCram ? detected.cram! : detected.bam!;
  const indexFile = isCram ? detected.crai! : detected.bai!;
  const formatLabel = isCram ? "CRAM" : "BAM";
  const steps = isCram ? 6 : 5; // +1 for variant screening

  // Step 1: Read index + header
  setProgress("Loading your data", `Loading ${formatLabel} index...`, 5, 1, steps);
  log(`Reading ${formatLabel} index...`);
  const indexBytes = new Uint8Array(await indexFile.arrayBuffer());
  log(`Index loaded: ${formatSize(indexBytes.length)}`);

  setProgress("Loading your data", `Loading ${formatLabel} header...`, 10, 1, steps);
  log(`Reading ${formatLabel} header...`);
  const headerSize = Math.min(dataFile.size, 4 * 1024 * 1024);
  const headerBytes = new Uint8Array(
    await dataFile.slice(0, headerSize).arrayBuffer(),
  );

  const regions: RegionInput[] = [
    { chrom: "chr17", start: 15229779, end: 15265326 },
    { chrom: "chr2", start: 50000000, end: 50500000 },
    { chrom: "chr17", start: 13500000, end: 21000000 },
    { chrom: "chr17", start: 17681458, end: 17811453 },
  ];

  // Step 2: Compute ranges
  setProgress("Preparing analysis", "Analysing index...", 15, 2, steps);
  log("Computing byte ranges from index...");
  const rangesReply = isCram
    ? wasm.required_ranges_cram(
        indexBytes,
        headerBytes,
        JSON.stringify(regions),
      )
    : wasm.required_ranges_bam(
        indexBytes,
        headerBytes,
        JSON.stringify(regions),
      );

  let reply: AnalyzeReply;
  let loadedData: Uint8Array | null = null; // for pileup reuse
  const config = JSON.stringify({
    chr2_control_start: 50000000,
    chr2_control_end: 50500000,
  });

  if (
    rangesReply.status === "ok" &&
    rangesReply.ranges &&
    rangesReply.ranges.length > 0
  ) {
    const MAX_CHUNK = 32 * 1024 * 1024;
    const rawRanges = rangesReply.ranges as { offset: number; length: number }[];
    const ranges: { offset: number; length: number }[] = [];
    for (const r of rawRanges) {
      let off = r.offset;
      let remaining = r.length;
      while (remaining > 0) {
        const chunk = Math.min(remaining, MAX_CHUNK);
        ranges.push({ offset: off, length: chunk });
        off += chunk;
        remaining -= chunk;
      }
    }

    const totalBytes = ranges.reduce((s, r) => s + r.length, 0);
    const pct = ((totalBytes / dataFile.size) * 100).toFixed(1);
    log(
      `Reading ${formatSize(totalBytes)} from ${formatSize(dataFile.size)} file (${pct}%)`,
    );

    // Step 3: Fetch byte ranges from the BAM/CRAM file.
    const t0 = performance.now();
    const combined = new Uint8Array(totalBytes);
    let pos = 0;

    for (let i = 0; i < ranges.length; i++) {
      const r = ranges[i]!;
      try {
        const buf = await sliceWithRetry(dataFile, r.offset, r.length);
        combined.set(new Uint8Array(buf), pos);
        pos += buf.byteLength;
        // Map fetch progress to 15-55% of overall bar
        const fetchPct = 15 + (pos / totalBytes) * 40;
        setProgress(
          "Loading genome data",
          `${formatSize(pos)} / ${formatSize(totalBytes)}`,
          fetchPct,
          3,
          steps,
        );
      } catch (e) {
        log(
          `FAILED at chunk ${i}: offset=${r.offset.toLocaleString()} len=${formatSize(r.length)}`,
        );
        log(`Error: ${e}`);
        throw e;
      }
    }
    const fetchMs = (performance.now() - t0).toFixed(0);
    log(`Data loaded in ${fetchMs} ms`);
    loadedData = combined;

    const rangesJson = JSON.stringify(
      ranges.map((r) => ({ file_offset: r.offset, length: r.length })),
    );

    // Step 4 (CRAM only): load reference FASTA, maps to 55-85%.
    if (isCram && detected.ref) {
      await loadReferenceFromFasta(detected.ref);
    }

    // SV analysis step (synchronous, blocks main thread).
    const svStep = isCram ? 5 : 4;
    setProgress(
      "Checking for duplications & deletions",
      "Scanning for large duplications & deletions (page may pause briefly)...",
      88,
      svStep,
      steps,
    );
    log("Running analysis...");
    await yieldToUI();

    const t1 = performance.now();
    reply = isCram
      ? wasm.analyze_cram_sparse(
          combined,
          rangesJson,
          indexBytes,
          cachedRefChr17!,
          cachedRefChr2!,
          config,
        )
      : wasm.analyze_bam_sparse(combined, rangesJson, indexBytes, config);
    const analyzeMs = (performance.now() - t1).toFixed(0);
    log(`Analysis complete in ${analyzeMs} ms`);
  } else {
    if (rangesReply.status === "error") {
      log(`Range query failed: ${rangesReply.error}. Reading full file...`);
    }
    setProgress("Loading genome data", "Loading full file...", 20, 3, steps);
    log(`Reading full ${formatLabel} file...`);
    const t0 = performance.now();
    const fullDataBytes = new Uint8Array(await dataFile.arrayBuffer());
    const fetchMs = (performance.now() - t0).toFixed(0);
    log(`File read in ${fetchMs} ms (${formatSize(fullDataBytes.length)})`);

    setProgress(
      "Checking for duplications & deletions",
      "Scanning for large duplications & deletions (page may pause briefly)...",
      88,
      isCram ? 5 : 4,
      steps,
    );
    log("Running analysis...");
    await yieldToUI();

    const t1 = performance.now();
    reply = wasm.analyze_bam(fullDataBytes, indexBytes, config);
    const analyzeMs = (performance.now() - t1).toFixed(0);
    log(`Analysis complete in ${analyzeMs} ms`);
    loadedData = fullDataBytes;
  }

  // Append summary from WASM progress events (skip walk_tick noise)
  if (reply.progress) {
    for (const evt of reply.progress) {
      const line = formatProgressEvent(evt);
      if (line) log(line);
    }
  }

  result.value = reply;

  // Detect biological sex from CRAI/BAI alignment spans for sex-aware
  // hemizygous calling on chrX/Y.
  let detectedSex = "unknown";
  if (isCram) {
    try {
      detectedSex = wasm.detect_sex_from_index(indexBytes, headerBytes);
      detectedSexValue.value = detectedSex;
      log(`Sex detection: ${detectedSex}`);
    } catch (e) {
      log(`Sex detection failed: ${e}`);
    }
  }

  // SNV pileup: process one chromosome at a time to keep memory bounded.
  // For each chrom: load reference from FASTA, load CRAM data, pileup, free.
  const snvStep = isCram ? 6 : 5;
  setProgress("Checking known disease-causing changes", "Preparing...", 92, snvStep, steps);
  log("Starting variant screening (chromosome by chromosome)...");

  try {
    const catalogRegionsJson: string = wasm.catalog_regions();
    let allRegions: RegionInput[] = JSON.parse(catalogRegionsJson);

    // Map ClinVar chromosome names to UCSC names used in CRAM/BAM
    const addPrefix = detectedBuild.value === "GRCh38" || detectedBuild.value === "GRCh37";
    if (addPrefix) {
      allRegions = allRegions.map((r) => {
        let chrom = r.chrom;
        if (!chrom.startsWith("chr")) chrom = `chr${chrom}`;
        // ClinVar uses "MT", UCSC/FASTA uses "chrM"
        if (chrom === "chrMT") chrom = "chrM";
        return { ...r, chrom };
      });
    }

    // Group regions by chromosome
    const byChrom = new Map<string, RegionInput[]>();
    for (const r of allRegions) {
      const list = byChrom.get(r.chrom) ?? [];
      list.push(r);
      byChrom.set(r.chrom, list);
    }

    // Natural chromosome order: 1,2,...,22,X,Y,MT
    const chromOrder = (c: string): number => {
      const n = c.replace("chr", "");
      if (n === "X") return 23;
      if (n === "Y") return 24;
      if (n === "MT" || n === "M" || n === "chrM") return 25;
      return parseInt(n, 10) || 99;
    };
    const chroms = [...byChrom.keys()].sort((a, b) => chromOrder(a) - chromOrder(b));
    log(`${allRegions.length} regions across ${chroms.length} chromosomes`);

    // Initialize per-chromosome screening state
    screeningProgress.chromosomes = chroms.map((c) => ({
      chrom: c,
      status: "pending" as ChromStatus,
      variantCount: 0,
      findings: 0,
    }));

    // Aggregate results across chromosomes
    const allFindings: CheckedVariant[] = [];
    let totalChecked = 0;
    let noCoverageCount = 0;
    let lowCoverageCount = 0;
    let chromsDone = 0;
    const refFile = detected.ref;
    const catalogMeta = JSON.parse(wasm.catalog_info());

    // Show initial pileupResult immediately so the UI renders
    pileupResult.value = {
      positive_findings: [],
      total_checked: 0,
      no_coverage_count: 0,
      low_coverage_count: 0,
      catalog_date: catalogMeta.generated_date ?? "",
      catalog_gene_count: catalogMeta.gene_count ?? 0,
      catalog_variant_count: catalogMeta.variant_count ?? 0,
    };

    for (const chrom of chroms) {
      chromsDone++;
      const chromRegions = byChrom.get(chrom)!;
      const pct = 92 + (chromsDone / chroms.length) * 7;
      screeningProgress.currentChrom = chrom;
      const chromState = screeningProgress.chromosomes.find((c) => c.chrom === chrom);
      if (chromState) chromState.status = "loading";
      setProgress(
        "Checking known disease-causing changes",
        `Chromosome ${chrom.replace("chr", "")} (${chromsDone}/${chroms.length})`,
        pct, snvStep, steps,
      );

      // 1. Load reference for this chromosome from FASTA
      let refSeq: Uint8Array | null = null;
      if (isCram && refFile) {
        // Strip "chr" prefix for FASTA scan (FASTA uses "chr1", "chr17")
        const fastaChromName = chrom;
        try {
          const chromPos = await findChromInFasta(refFile, fastaChromName);
          refSeq = await readChromSequence(refFile, fastaChromName, chromPos.seqStart, chromPos.seqEnd);
        } catch {
          log(`  ${chrom}: not in reference FASTA, skipping`);
          noCoverageCount += chromRegions.length;
          if (chromState) chromState.status = "skipped";
          continue;
        }
      }

      // 2. Compute CRAM byte ranges for this chromosome's regions (tight margin)
      const rangesReply = isCram
        ? wasm.required_ranges_cram_tight(indexBytes, headerBytes, JSON.stringify(chromRegions))
        : wasm.required_ranges_bam(indexBytes, headerBytes, JSON.stringify(chromRegions));

      if (rangesReply.status !== "ok" || !rangesReply.ranges || rangesReply.ranges.length === 0) {
        log(`  ${chrom}: no data in file, skipping`);
        if (chromState) chromState.status = "skipped";
        continue;
      }

      // 3. File.slice() the CRAM data for this chromosome
      const MAX_CHUNK = 32 * 1024 * 1024;
      const rawRanges = rangesReply.ranges as { offset: number; length: number }[];
      const chromRanges: { offset: number; length: number }[] = [];
      for (const r of rawRanges) {
        let off = r.offset;
        let remaining = r.length;
        while (remaining > 0) {
          const chunk = Math.min(remaining, MAX_CHUNK);
          chromRanges.push({ offset: off, length: chunk });
          off += chunk;
          remaining -= chunk;
        }
      }

      const chromTotal = chromRanges.reduce((s, r) => s + r.length, 0);
      const chromData = new Uint8Array(chromTotal);
      let pos = 0;
      for (const r of chromRanges) {
        const buf = await sliceWithRetry(dataFile, r.offset, r.length);
        chromData.set(new Uint8Array(buf), pos);
        pos += buf.byteLength;
      }

      // 4. Run pileup for this chromosome
      const chromRangesJson = JSON.stringify(
        chromRanges.map((r) => ({ file_offset: r.offset, length: r.length })),
      );

      try {
        const chromResult = isCram
          ? wasm.check_variants_cram_chrom(chromData, chromRangesJson, indexBytes, refSeq!, chrom, detectedSex)
          : wasm.check_variants_bam_sparse(chromData, chromRangesJson, indexBytes);

        if (chromResult?.positive_findings) {
          allFindings.push(...chromResult.positive_findings);
        }
        totalChecked += chromResult?.total_checked ?? 0;
        noCoverageCount += chromResult?.no_coverage_count ?? 0;
        lowCoverageCount += chromResult?.low_coverage_count ?? 0;

        const found = chromResult?.positive_findings?.length ?? 0;
        const checked = chromResult?.total_checked ?? 0;
        if (chromState) {
          chromState.status = "done";
          chromState.variantCount = checked;
          chromState.findings = found;
        }
        if (found > 0) {
          log(`  ${chrom}: ${found} findings from ${checked} variants`);
        }
      } catch (e) {
        log(`  ${chrom}: error: ${e}`);
        if (chromState) chromState.status = "error";
      }

      // Update pileupResult progressively so the UI refreshes live
      pileupResult.value = {
        positive_findings: allFindings,
        total_checked: totalChecked,
        no_coverage_count: noCoverageCount,
        low_coverage_count: lowCoverageCount,
        catalog_date: catalogMeta.generated_date ?? "",
        catalog_gene_count: catalogMeta.gene_count ?? 0,
        catalog_variant_count: catalogMeta.variant_count ?? 0,
      };

      await yieldToUI();
    }

    screeningProgress.currentChrom = "";
    log(`Variant screening complete: ${allFindings.length} findings from ${totalChecked} variants across ${chroms.length} chromosomes`);
  } catch (e) {
    log(`Variant screening error: ${e}`);
  }

  setProgress("Complete", "", 100);
}

function formatProgressEvent(evt: ProgressEvent): string | null {
  switch (evt.kind) {
    case "walk_tick":
      return null; // suppress noisy tick events
    case "header_verified":
      return "Header verified";
    case "walk_start":
      return `Scanning ${evt.phase ?? ""}`;
    case "walk_done":
      return `${evt.phase ?? ""}: ${(evt.records_seen ?? 0).toLocaleString()} records`;
    case "depth_ratio_done":
      return `Depth ratio: ${evt.ratio?.toFixed(3) ?? "?"} (CN=${evt.estimated_cn ?? "?"})`;
    case "boundary_scan_done":
      return `Boundary scan: ${evt.windows ?? 0} windows`;
    case "interpret_done":
      return `Result: ${evt.confidence ?? "?"}, subtype=${evt.subtype_call ?? "none"}`;
    default:
      return evt.kind;
  }
}

export function useAnalysis() {
  async function analyze(wasm: WasmModule, files: DetectedFiles) {
    state.value = "loading";
    logs.value = [];
    result.value = null;
    pileupResult.value = null;
    setProgress("Starting", "", 0, 0, 0);
    try {
      state.value = "analyzing";
      await runAnalysis(wasm, files);
      state.value = "done";
    } catch (e) {
      log(`ERROR: ${e}`);
      state.value = "error";
    }
  }

  function reset() {
    state.value = "idle";
    logs.value = [];
    result.value = null;
    pileupResult.value = null;
    screeningProgress.chromosomes = [];
    screeningProgress.currentChrom = "";
    setProgress("", "", 0, 0, 0);
    cachedRefChr17 = null;
    cachedRefChr2 = null;
  }

  async function detectFileBuild(
    wasm: WasmModule,
    file: File,
    isCram: boolean,
  ): Promise<string | null> {
    try {
      const headerSize = Math.min(file.size, 4 * 1024 * 1024);
      const headerBuf = new Uint8Array(
        await file.slice(0, headerSize).arrayBuffer(),
      );
      const build = wasm.detect_build(headerBuf, isCram);
      detectedBuild.value = build;
      return build;
    } catch {
      return null;
    }
  }

  const isRunning = computed(
    () => state.value === "loading" || state.value === "analyzing",
  );

  return {
    state,
    logs,
    result,
    pileupResult,
    screeningProgress,
    progress,
    detectedBuild,
    detectedSex: detectedSexValue,
    isRunning,
    analyze,
    reset,
    detectFileBuild,
  };
}
