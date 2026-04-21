<script setup lang="ts">
import { computed } from "vue";
import type {
  AnalyzeReply,
  CheckedVariant,
  PileupResult,
  WasmModule,
  ScreeningManifest,
  SvScreeningItem,
  SnvScreeningItem,
} from "../types/wasm";
import type { ScreeningProgress } from "../composables/useAnalysis";
import { formatSize } from "../composables/useAnalysis";
import { useWasm } from "../composables/useWasm";
import ChromosomeDiagram from "./ChromosomeDiagram.vue";
import SnvFindingCard from "./SnvFindingCard.vue";

const props = defineProps<{
  svResult: AnalyzeReply;
  pileupResult: PileupResult | null;
  screening: ScreeningProgress;
  detectedBuild: string | null;
  detectedSex: string;
  wasm: WasmModule | null;
}>();

const { getReferenceData } = useWasm();

const interp = computed(() => props.svResult.interpretation);
const depth = computed(() => props.svResult.depth);
const boundary = computed(() => props.svResult.boundary);
const runs = computed(() => boundary.value?.runs ?? []);
const run = computed(() => runs.value[0] ?? null);

const refData = computed(() =>
  getReferenceData(props.detectedBuild ?? "GRCh38"),
);
const chr17Genes = computed(() =>
  (refData.value?.genes ?? []).filter((g) => g.region.chrom === "chr17"),
);
const chr17Reps = computed(() => refData.value?.reps ?? []);

// Manifest from Rust - the single source of truth
const manifest = computed<ScreeningManifest | null>(() => {
  if (!props.wasm) return null;
  try {
    return JSON.parse(
      props.wasm.screening_manifest(props.detectedBuild ?? "GRCh38"),
    );
  } catch {
    return null;
  }
});

// SV results mapped to manifest items
const svResults = computed(() => {
  if (!manifest.value || !interp.value) return [];
  const subtype = interp.value.subtype_call?.toUpperCase() ?? null;
  return manifest.value.sv_screening.map((item) => ({
    ...item,
    detected: subtype === item.short_name.toUpperCase(),
  }));
});

// SNV positive findings. A heterozygous variant in a recessive gene is a
// carrier-status finding; everything else is disease-causing.
function isCarrierFinding(f: CheckedVariant): boolean {
  const isHet = typeof f.call === "object" && "Heterozygous" in f.call;
  return isHet && (f.inheritance === "AutosomalRecessive" || f.inheritance === "Both");
}

const snvFindings = computed<CheckedVariant[]>(
  () => props.pileupResult?.positive_findings ?? [],
);
const diseaseFindings = computed(() => snvFindings.value.filter((f) => !isCarrierFinding(f)));
const carrierFindings = computed(() => snvFindings.value.filter((f) => isCarrierFinding(f)));
const isScreening = computed(() => props.screening.currentChrom !== "");

// Per-gene finding severity: "disease" if any non-carrier, else "carrier".
const geneFindingSeverity = computed(() => {
  const map = new Map<string, "disease" | "carrier">();
  for (const f of snvFindings.value) {
    if (!isCarrierFinding(f)) {
      map.set(f.gene, "disease");
    } else if (!map.has(f.gene)) {
      map.set(f.gene, "carrier");
    }
  }
  return map;
});
const chromFindingSeverity = computed(() => {
  const map = new Map<string, "disease" | "carrier">();
  for (const f of snvFindings.value) {
    const chrom = `chr${f.chrom}`;
    if (!isCarrierFinding(f)) {
      map.set(chrom, "disease");
    } else if (!map.has(chrom)) {
      map.set(chrom, "carrier");
    }
  }
  return map;
});

// Chromosomes that have been processed (for coloring gene chips)
const processedChroms = computed(() => {
  const set = new Set<string>();
  for (const c of props.screening.chromosomes) {
    if (c.status === "done" || c.status === "skipped") {
      set.add(c.chrom);
    }
  }
  return set;
});

function geneChipClass(gene: { gene: string; chrom: string }): string {
  const severity = geneFindingSeverity.value.get(gene.gene);
  if (severity === "disease") return "gene-finding";
  if (severity === "carrier") return "gene-carrier";
  // During screening: pulse genes on the current chromosome
  if (isScreening.value) {
    const currentChrom = props.screening.currentChrom;
    // Match "17" to "chr17" etc.
    const geneChrom = gene.chrom.startsWith("chr") ? gene.chrom : `chr${gene.chrom}`;
    if (geneChrom === currentChrom) return "gene-scanning";
    // Already-processed chromosomes: check if done
    if (processedChroms.value.has(geneChrom)) return "gene-clear";
    return "";
  }
  if (props.pileupResult) return "gene-clear";
  return "";
}

// Verdict
const verdictClass = computed(() => {
  const conf = interp.value?.confidence;
  const sub = interp.value?.subtype_call;
  const hasSnvFindings = snvFindings.value.length > 0;
  if (conf === "Refused") return "verdict-refused";
  if (sub || hasSnvFindings) return "verdict-finding";
  if (interp.value?.copy_number !== null && interp.value?.copy_number !== 2)
    return "verdict-atypical";
  if (isScreening.value) return "verdict-atypical"; // amber while still scanning
  return "verdict-normal";
});

const verdictHeadline = computed(() => {
  const sub = interp.value?.subtype_call;
  const snvCount = snvFindings.value.length;
  if (interp.value?.confidence === "Refused") return "Insufficient Coverage";
  if (sub && snvCount > 0) return `Findings: ${sub} + ${snvCount} gene change(s)`;
  if (sub) return `Finding: ${sub}`;
  if (snvCount > 0) return `${snvCount} Disease-Causing Change(s) Found`;
  if (interp.value?.copy_number !== null && interp.value?.copy_number !== 2)
    return "Unusual Finding";
  return isScreening.value ? "Duplication/deletion check done - gene screening in progress..." : "No Findings Detected";
});

// Next steps
const nextSteps = computed(() => {
  const sub = interp.value?.subtype_call;
  const snvCount = snvFindings.value.length;
  if (interp.value?.confidence === "Refused")
    return "Insufficient coverage for reliable screening. At least 15x whole-genome sequencing is needed.";
  if (sub || snvCount > 0)
    return "This screening result suggests possible findings. This is NOT a diagnosis. Please take this result to a clinical geneticist or neurologist for confirmation.";
  if (isScreening.value)
    return "Duplication/deletion screening found no findings. Gene-level screening is still in progress...";
  return "No disease-causing changes detected in the screened conditions. This does not completely rule out hereditary neuropathy - changes not yet in medical databases or in genes outside our panel are not detected.";
});

function downloadReport() {
  if (!interp.value || !depth.value) return;
  const build = props.detectedBuild ?? "unknown";
  const lines = [
    "Neuropathy DNA Scanner screening report",
    "=".repeat(40),
    "",
    "NOT a medical diagnosis. Confirm with clinical genetic testing.",
    "",
    `Reference build: ${build}`,
    `Copy number: ${interp.value.copy_number ?? "not determined"}`,
    `Confidence: ${interp.value.confidence}`,
    `SV subtype: ${interp.value.subtype_call ?? "none"}`,
    `SNV findings: ${snvFindings.value.length}`,
    "",
    `Depth ratio: ${depth.value.ratio.toFixed(3)}`,
    `  PMP22: ${depth.value.pmp22_mean_depth.toFixed(1)}x`,
    `  Autosomal: ${depth.value.autosomal_mean_depth.toFixed(1)}x`,
    "",
    interp.value.plain_language,
  ];
  if (snvFindings.value.length > 0) {
    lines.push("", "SNV findings:");
    for (const v of snvFindings.value) {
      lines.push(`  ${v.gene} ${v.chrom}:${v.pos} ${v.ref_allele}>${v.alt_allele} (ClinVar:${v.clinvar_id})`);
    }
  }
  lines.push("", "--- Disclaimer ---", "NOT a medical device. No liability assumed.");
  const blob = new Blob([lines.join("\n")], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "cmt_result.txt";
  a.click();
  URL.revokeObjectURL(url);
}
</script>

<template>
  <section class="results">
    <!-- Reminder -->
    <div class="result-reminder">
      This is a screening tool, not a diagnosis. Discuss any result with a
      medical professional.
    </div>

    <!-- Verdict -->
    <div class="verdict" :class="verdictClass">
      <h2>{{ verdictHeadline }}</h2>
      <p>{{ nextSteps }}</p>
    </div>

    <!-- ============================================ -->
    <!-- SCREENING CHECKLIST - all conditions tested  -->
    <!-- ============================================ -->
    <div v-if="manifest" class="checklist">
      <!-- SV: large duplications and deletions -->
      <h3>Large Duplications and Deletions</h3>
      <p class="section-note">
        Checks whether parts of your DNA are duplicated or missing.
        These are the most common causes of inherited neuropathies.
      </p>
      <div class="checklist-items">
        <div
          v-for="s in svResults"
          :key="s.short_name"
          class="check-row"
          :class="{ detected: s.detected }"
        >
          <span class="check-icon">{{ s.detected ? "FOUND" : "clear" }}</span>
          <span class="check-label">
            <strong>{{ s.short_name }}</strong> - {{ s.name }}
            <span v-if="s.extra_check" class="bonus-tag">extra</span>
          </span>
          <a
            v-if="s.omim_url"
            :href="s.omim_url"
            target="_blank"
            rel="noopener"
            class="check-link"
            >OMIM</a
          >
          <a
            v-if="s.genereviews_url"
            :href="s.genereviews_url"
            target="_blank"
            rel="noopener"
            class="check-link"
            >GeneReviews</a
          >
        </div>
      </div>

      <!-- ============================================ -->
      <!-- POINT VARIANT SCREENING                      -->
      <!-- ============================================ -->
      <h3>
        Known Disease-Causing DNA Changes
        <span v-if="isScreening" class="screening-spinner"></span>
      </h3>
      <p class="section-note">
        Checks {{ manifest.catalog_gene_count }} genes for
        {{ manifest.catalog_variant_count.toLocaleString() }} known disease-causing
        changes, one chromosome at a time.
      </p>

      <!-- SNV positive findings: disease-causing first, then carriers -->
      <template v-if="diseaseFindings.length > 0">
        <h4 class="findings-heading">Disease-Causing Changes Detected</h4>
        <SnvFindingCard
          v-for="v in diseaseFindings"
          :key="`${v.chrom}:${v.pos}:${v.alt_allele}`"
          :variant="v"
        />
      </template>
      <template v-if="carrierFindings.length > 0">
        <h4 class="findings-heading carrier-heading">Carrier Status</h4>
        <p class="section-note">
          You carry one copy of a recessive variant. This does not cause disease
          on its own but may be relevant for family planning.
        </p>
        <SnvFindingCard
          v-for="v in carrierFindings"
          :key="`${v.chrom}:${v.pos}:${v.alt_allele}`"
          :variant="v"
        />
      </template>

      <div v-if="pileupResult && !isScreening && snvFindings.length === 0" class="snv-clear">
        No known disease-causing DNA changes detected
      </div>

      <!-- Chromosome progress -->
      <div v-if="screening.chromosomes.length > 0" class="chrom-section">
        <p class="chrom-label">Chromosomes checked:</p>
        <div class="chrom-grid">
          <div
            v-for="c in screening.chromosomes"
            :key="c.chrom"
            class="chrom-chip"
            :class="[c.status, chromFindingSeverity.get(c.chrom) ?? '']"
            :title="c.status === 'skipped'
              ? `Chromosome ${c.chrom.replace('chr', '')}: skipped (not in reference)`
              : `Chromosome ${c.chrom.replace('chr', '')}: ${c.variantCount} changes checked, ${c.findings} found`"
          >
            <span>{{ c.chrom.replace("chr", "") }}</span>
            <span v-if="c.status === 'loading'" class="chip-spin"></span>
          </div>
        </div>
      </div>

      <p v-if="detectedSex !== 'unknown'" class="sex-note">
        Detected sex: {{ detectedSex }}. chrX variants assessed with
        {{ detectedSex === 'male' ? 'hemizygous' : 'diploid' }} thresholds.
      </p>

      <!-- CMT core genes -->
      <details v-if="manifest.cmt_core_genes.length" class="gene-list-details" open>
        <summary>
          CMT / Neuropathy Genes
          ({{ manifest.cmt_core_genes.length }} genes)
        </summary>
        <div class="gene-grid">
          <span
            v-for="g in manifest.cmt_core_genes"
            :key="g.gene"
            :class="['gene-chip', 'cmt-core', geneChipClass(g)]"
            :title="`${g.gene}: ${g.variant_count} known disease-causing changes in catalog, 0 detected in this sample. Condition: ${g.condition}. Inheritance: ${g.inheritance}.`"
          >
            {{ g.gene }}
            <span class="gene-chip-count">{{ g.variant_count }}</span>
          </span>
        </div>
      </details>

      <!-- Extended panel genes -->
      <details v-if="manifest.extended_panel_genes.length" class="gene-list-details">
        <summary>
          Other Neuropathy-Related Genes
          ({{ manifest.extended_panel_genes.length }} genes)
        </summary>
        <p class="panel-note">
          Genes linked to conditions that can affect the peripheral nerves.
          Includes rare metabolic, degenerative, and syndromic conditions.
        </p>
        <div class="gene-grid">
          <span
            v-for="g in manifest.extended_panel_genes"
            :key="g.gene"
            :class="['gene-chip', geneChipClass(g)]"
            :title="`${g.gene}: ${g.variant_count} known disease-causing changes in catalog, 0 detected. Condition: ${g.condition}. Inheritance: ${g.inheritance}.`"
          >
            {{ g.gene }}
            <span class="gene-chip-count">{{ g.variant_count }}</span>
          </span>
        </div>
      </details>

      <!-- Variant count + filter explanation -->
      <div v-if="pileupResult" class="catalog-stats">
        <span>
          {{ manifest.catalog_variant_count.toLocaleString() }} known changes in catalog
        </span>
        <span>
          {{ manifest.catalog_snv_count.toLocaleString() }} single-letter changes checked
        </span>
        <span v-if="manifest.catalog_skipped_indels > 0" class="stat-note">
          {{ manifest.catalog_skipped_indels.toLocaleString() }} insertions/deletions not yet supported
        </span>
        <span v-if="pileupResult.no_coverage_count > 0" class="stat-note">
          {{ pileupResult.no_coverage_count.toLocaleString() }} positions had no data
        </span>
      </div>

      <p class="catalog-meta">
        Catalog: ClinVar {{ manifest.catalog_date }}
        ({{ manifest.catalog_variant_count.toLocaleString() }} variants,
        {{ manifest.catalog_gene_count }} genes).
        Filter: {{ manifest.catalog_filter }}.
        Source: {{ manifest.catalog_source }}.
      </p>

      <!-- Narrative linking SV + SNV results -->
      <div v-if="pileupResult && !isScreening" class="narrative-link">
        <template v-if="interp?.subtype_call && snvFindings.length === 0">
          {{ interp.subtype_call }} was detected via copy-number analysis.
          Separately, targeted genotyping found no additional pathogenic
          point variants in {{ manifest.catalog_gene_count }} neuropathy-associated genes.
          This rules out combined causes but does not change the
          {{ interp.subtype_call }} finding.
        </template>
        <template v-else-if="!interp?.subtype_call && snvFindings.length === 0">
          No structural variants or known pathogenic point variants were detected.
          This does not exclude hereditary neuropathy - novel variants and
          conditions outside the catalog are not covered.
        </template>
      </div>

      <p v-if="pileupResult" class="catalog-meta">
        Catalog: ClinVar {{ pileupResult.catalog_date }}
      </p>
    </div>

    <!-- ============================================ -->
    <!-- TECHNICAL DETAILS (expandable)               -->
    <!-- ============================================ -->
    <details class="expandable">
      <summary>Technical details</summary>
      <table class="detail-table">
        <tr>
          <td>Copy number</td>
          <td><strong>{{ interp?.copy_number ?? "--" }}</strong> (normal = 2)</td>
        </tr>
        <tr>
          <td>Depth ratio</td>
          <td>
            {{ depth?.ratio.toFixed(3) }} (PMP22: {{ depth?.pmp22_mean_depth.toFixed(1) }}x,
            autosomal: {{ depth?.autosomal_mean_depth.toFixed(1) }}x)
          </td>
        </tr>
        <tr>
          <td>Confidence</td>
          <td>{{ interp?.confidence === "Full" ? "High" : interp?.confidence }}</td>
        </tr>
        <tr v-if="run">
          <td>Detected region</td>
          <td>
            chr17:{{ run.start.toLocaleString() }}-{{ run.end.toLocaleString() }}
            ({{ formatSize(run.length) }}, {{ run.direction }})
          </td>
        </tr>
        <tr v-if="run">
          <td>Signal quality</td>
          <td>{{ run.bridged_gaps === 0 ? "Clean" : `${run.bridged_gaps} gap(s) bridged` }}</td>
        </tr>
        <tr v-if="depth?.rai1_mean_depth">
          <td>RAI1 depth</td>
          <td>{{ depth.rai1_mean_depth.toFixed(1) }}x (CN={{ depth.rai1_estimated_cn }})</td>
        </tr>
        <tr><td>Build</td><td>{{ detectedBuild ?? "unknown" }}</td></tr>
        <tr><td>Analysis time</td><td>{{ svResult.elapsed_ms ?? "?" }} ms</td></tr>
      </table>

      <template v-if="run && chr17Reps.length">
        <h4 class="sub-heading">Position on chromosome 17</h4>
        <ChromosomeDiagram :reps="chr17Reps" :genes="chr17Genes" :run="run" />
      </template>
    </details>

    <!-- Download -->
    <div class="download-row">
      <button class="btn-download" @click="downloadReport">
        Download report
      </button>
    </div>

    <div class="result-footer">{{ interp?.plain_language }}</div>
  </section>
</template>

<style scoped>
.results {
  margin: 1rem 0;
}

.result-reminder {
  background: var(--amber-light);
  border: 1px solid var(--amber);
  border-radius: 4px;
  padding: 0.5rem 0.75rem;
  font-size: 0.85rem;
  margin-bottom: 1rem;
}

/* Verdict */
.verdict { padding: 1rem; border-radius: 6px; margin-bottom: 1rem; }
.verdict h2 { margin: 0 0 0.5rem; }
.verdict p { font-size: 0.9rem; }
.verdict-normal { background: var(--green-light); border-left: 4px solid var(--green); }
.verdict-normal h2 { color: var(--green); }
.verdict-finding { background: var(--red-light); border-left: 4px solid var(--red); }
.verdict-finding h2 { color: var(--red); }
.verdict-atypical { background: var(--amber-light); border-left: 4px solid var(--amber); }
.verdict-atypical h2 { color: var(--amber); }
.verdict-refused { background: var(--code-bg); border-left: 4px solid var(--text-light); }

/* Checklist */
.checklist { margin: 1rem 0; }
.checklist h3 {
  font-size: 1rem;
  margin: 1rem 0 0.5rem;
  display: flex;
  align-items: center;
  gap: 0.5rem;
}
.section-note { font-size: 0.8rem; color: var(--text-light); margin: -0.3rem 0 0.5rem; }
.chrom-section { margin: 0.5rem 0; }
.chrom-label { font-size: 0.75rem; color: var(--text-light); margin-bottom: 0.25rem; }
.sex-note { font-size: 0.72rem; color: var(--text-light); margin: 0.25rem 0 0.5rem; }
.checklist-items { display: flex; flex-direction: column; gap: 3px; margin-bottom: 0.75rem; }
.check-row {
  display: flex;
  align-items: center;
  gap: 0.5rem;
  padding: 0.3rem 0.5rem;
  border-radius: 4px;
  font-size: 0.85rem;
  background: var(--green-light);
}
.check-row.detected { background: var(--red-light); }
.check-icon {
  font-size: 0.6rem;
  font-weight: 700;
  min-width: 40px;
  text-align: center;
}
.check-row:not(.detected) .check-icon { color: var(--green); }
.check-row.detected .check-icon { color: white; background: var(--red); border-radius: 3px; padding: 1px 4px; }
.check-label { flex: 1; }
.bonus-tag {
  font-size: 0.6rem;
  font-weight: normal;
  color: var(--text-light);
  background: var(--code-bg);
  border: 1px solid var(--border);
  border-radius: 3px;
  padding: 0 4px;
  margin-left: 0.4rem;
  vertical-align: middle;
}
.check-link {
  font-size: 0.65rem;
  color: var(--accent);
  border: 1px solid var(--accent);
  border-radius: 3px;
  padding: 0 4px;
  text-decoration: none;
  opacity: 0.7;
}
.check-link:hover { opacity: 1; background: #eff6ff; }

/* Chromosome chips */
.chrom-grid { display: flex; flex-wrap: wrap; gap: 3px; margin: 0.5rem 0; }
.chrom-chip {
  display: inline-flex; align-items: center; gap: 2px;
  padding: 2px 5px; border-radius: 3px;
  font-size: 0.65rem; font-weight: 600;
  border: 1px solid var(--border); background: var(--code-bg); color: var(--text-light);
}
.chrom-chip.done { background: var(--green-light); border-color: var(--green); color: var(--green); }
.chrom-chip.disease { background: #fef2f2; border-color: #dc2626; color: #dc2626; }
.chrom-chip.carrier { background: #fffbeb; border-color: #d97706; color: #d97706; }
.chrom-chip.loading {
  background: #eff6ff; border-color: var(--accent); color: var(--accent);
  animation: pulse-gene 1.2s ease-in-out infinite;
}
.chrom-chip.skipped { opacity: 0.4; }
.chip-spin {
  width: 7px; height: 7px;
  border: 1.5px solid var(--border); border-top-color: var(--accent);
  border-radius: 50%; animation: spin 0.8s linear infinite;
}
.screening-spinner {
  width: 12px; height: 12px;
  border: 2px solid var(--border); border-top-color: var(--accent);
  border-radius: 50%; animation: spin 0.8s linear infinite;
  display: inline-block;
}
@keyframes spin { to { transform: rotate(360deg); } }

/* Gene list */
.gene-list-details { margin: 0.5rem 0; }
.gene-list-details summary { font-size: 0.8rem; cursor: pointer; color: var(--text-light); }
.gene-grid { display: flex; flex-wrap: wrap; gap: 3px; margin-top: 0.5rem; }
.gene-chip {
  font-size: 0.65rem; padding: 1px 5px; border-radius: 3px;
  background: var(--code-bg); border: 1px solid var(--border);
  font-style: italic;
}
.gene-chip-count {
  font-size: 0.55rem; color: var(--text-light); font-style: normal;
  margin-left: 2px;
}

/* SNV findings */
.findings-heading { font-size: 0.95rem; margin: 1rem 0 0.5rem; color: var(--red); }
.findings-heading.carrier-heading { color: #d97706; }
.snv-clear {
  background: var(--green-light); border: 1px solid var(--green);
  border-radius: 4px; padding: 0.5rem 0.75rem;
  font-size: 0.85rem; color: var(--green); font-weight: 600;
  margin: 0.5rem 0;
}
.gene-chip.cmt-core {
  background: #eff6ff;
  border-color: var(--accent);
}
.gene-chip.gene-clear {
  background: var(--green-light);
  border-color: var(--green);
  color: var(--green);
}
.gene-chip.gene-finding {
  background: var(--red-light);
  border-color: var(--red);
  color: var(--red);
  font-weight: 700;
}
.gene-chip.gene-carrier {
  background: #fffbeb;
  border-color: #d97706;
  color: #d97706;
  font-weight: 700;
}
.gene-chip.gene-scanning {
  background: #eff6ff;
  border-color: var(--accent);
  color: var(--accent);
  animation: pulse-gene 1.2s ease-in-out infinite;
}
@keyframes pulse-gene {
  0%, 100% { opacity: 0.5; transform: scale(1); }
  50% { opacity: 1; transform: scale(1.05); }
}
.panel-note {
  font-size: 0.75rem; color: var(--text-light); margin: 0.3rem 0;
}
.catalog-stats {
  font-size: 0.75rem; color: var(--text-light); margin: 0.5rem 0;
  display: flex; flex-wrap: wrap; gap: 0.5rem;
}
.stat-note { opacity: 0.8; }
.catalog-meta {
  font-size: 0.7rem; color: var(--text-light); margin-top: 0.3rem;
  line-height: 1.4;
}
.narrative-link {
  font-size: 0.85rem; color: var(--text-light);
  margin: 1rem 0; padding: 0.5rem 0.75rem;
  background: var(--code-bg); border-radius: 4px;
  line-height: 1.5;
}

/* Expandable */
.expandable { margin: 1rem 0; }
.expandable summary { cursor: pointer; font-weight: 600; padding: 0.3rem 0; }
.sub-heading { margin: 1rem 0 0.3rem; font-size: 0.95rem; }

.download-row { margin: 1rem 0; }
.result-footer {
  margin-top: 1rem; padding-top: 0.75rem;
  border-top: 1px solid var(--border);
  font-size: 0.85rem; color: var(--text-light);
}
</style>
