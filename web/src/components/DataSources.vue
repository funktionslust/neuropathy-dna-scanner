<script setup lang="ts">
import { computed } from "vue";
import { useWasm } from "../composables/useWasm";
import type { WasmModule, ReferenceDataset, ScreeningManifest } from "../types/wasm";

const props = defineProps<{
  wasm: WasmModule | null;
  detectedBuild: string | null;
}>();

const { getReferenceData } = useWasm();

const data = computed<ReferenceDataset | null>(() => {
  if (!props.wasm) return null;
  return getReferenceData(props.detectedBuild ?? "GRCh38");
});

const manifest = computed<ScreeningManifest | null>(() => {
  if (!props.wasm) return null;
  try {
    return JSON.parse(props.wasm.screening_manifest(props.detectedBuild ?? "GRCh38"));
  } catch { return null; }
});

function downloadCatalog() {
  if (!props.wasm || !manifest.value) return;
  // Get the full variant-level CSV from WASM (22,844 rows)
  const csv: string = props.wasm.catalog_csv();
  const blob = new Blob([csv], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `neuropathy_variant_catalog_${manifest.value.catalog_date}.csv`;
  a.click();
  URL.revokeObjectURL(url);
}
</script>

<template>
  <section v-if="data" class="data-sources">
    <details>
      <summary><h2 class="inline-heading">Data sources and coordinates</h2></summary>

      <p class="note">
        Build: {{ data.build }} | Source:
        <a :href="data.build_source" target="_blank" rel="noopener">{{
          data.build_source
        }}</a>
      </p>

      <!-- Genes -->
      <h4>Gene positions</h4>
      <table class="detail-table">
        <tr>
          <th>Gene</th>
          <th>Region</th>
          <th>Source</th>
        </tr>
        <tr v-for="g in data.genes" :key="g.symbol">
          <td>
            <strong>{{ g.symbol }}</strong>  - {{ g.name }}
          </td>
          <td>
            {{ g.region.chrom }}:{{
              g.region.start.value.toLocaleString()
            }}-{{ g.region.end.value.toLocaleString() }}
            <a
              :href="g.region.start.source_url"
              target="_blank"
              rel="noopener"
              class="src-badge"
              title="Verify this coordinate"
            >src</a>
          </td>
          <td>
            <a
              :href="g.region.start.source_url"
              target="_blank"
              rel="noopener"
              >{{ g.region.start.source }}</a
            >
          </td>
        </tr>
      </table>

      <!-- REPs -->
      <h4>Segmental duplication elements (REPs)</h4>
      <table class="detail-table">
        <tr>
          <th>Element</th>
          <th>Region</th>
          <th>Size</th>
          <th>Publication</th>
        </tr>
        <tr v-for="r in data.reps" :key="r.name">
          <td>{{ r.name }}</td>
          <td>
            {{ r.region.chrom }}:{{
              r.region.start.value.toLocaleString()
            }}-{{ r.region.end.value.toLocaleString() }}
            <a
              :href="r.publication_url"
              target="_blank"
              rel="noopener"
              class="src-badge"
              title="Verify this coordinate"
            >src</a>
          </td>
          <td>{{ r.size_bp.toLocaleString() }} bp</td>
          <td>
            <a :href="r.publication_url" target="_blank" rel="noopener">{{
              r.publication
            }}</a>
          </td>
        </tr>
      </table>

      <!-- SV Syndromes (from reference_data.rs) -->
      <h4>Large duplications and deletions</h4>
      <table class="detail-table">
        <tr>
          <th>Condition</th>
          <th>Mechanism</th>
          <th>Expected CN</th>
          <th>OMIM</th>
          <th>ICD-10</th>
        </tr>
        <tr v-for="s in data.syndromes" :key="s.short_name">
          <td>
            <strong>{{ s.short_name }}</strong> - {{ s.name }}
          </td>
          <td>{{ s.mechanism }}</td>
          <td>{{ s.expected_cn }}</td>
          <td>
            <a :href="s.omim_url" target="_blank" rel="noopener"
              >#{{ s.omim }}</a
            >
          </td>
          <td>{{ s.icd10 || "--" }}</td>
        </tr>
      </table>

      <!-- SNV Genes (from generated ClinVar catalog) -->
      <template v-if="manifest">
        <h4>
          Known disease-causing gene changes
          ({{ manifest.catalog_variant_count.toLocaleString() }} variants,
          {{ manifest.catalog_gene_count }} genes)
        </h4>
        <p class="note">
          Source: {{ manifest.catalog_source }}.
          Filter: {{ manifest.catalog_filter }}.
          Catalog date: {{ manifest.catalog_date }}.
        </p>

        <h5>CMT / Neuropathy core genes ({{ manifest.cmt_core_genes.length }})</h5>
        <table class="detail-table">
          <tr>
            <th>Gene</th>
            <th>Chromosome</th>
            <th>Condition</th>
            <th>Inheritance</th>
            <th>Variants</th>
          </tr>
          <tr v-for="g in manifest.cmt_core_genes" :key="g.gene">
            <td >
              <strong><em>{{ g.gene }}</em></strong>
              <a
                :href="`https://www.ncbi.nlm.nih.gov/clinvar/?term=${g.gene}[gene]+pathogenic[clinsig]`"
                target="_blank" rel="noopener" class="src-badge src-right"
                title="View pathogenic variants for this gene in ClinVar"
              >ClinVar</a>
            </td>
            <td>{{ g.chrom }}</td>
            <td>{{ g.condition }}</td>
            <td>{{ g.inheritance }}</td>
            <td>{{ g.variant_count }}</td>
          </tr>
        </table>

        <details>
          <summary class="note">
            Other neuropathy-related genes ({{ manifest.extended_panel_genes.length }})
          </summary>
          <table class="detail-table">
            <tr>
              <th>Gene</th>
              <th>Chr</th>
              <th>Condition</th>
              <th>Inheritance</th>
              <th>Variants</th>
            </tr>
            <tr v-for="g in manifest.extended_panel_genes" :key="g.gene">
              <td >
                <em>{{ g.gene }}</em>
                <a
                  :href="`https://www.ncbi.nlm.nih.gov/clinvar/?term=${g.gene}[gene]+pathogenic[clinsig]`"
                  target="_blank" rel="noopener" class="src-badge src-right"
                >ClinVar</a>
              </td>
              <td>{{ g.chrom }}</td>
              <td>{{ g.condition }}</td>
              <td>{{ g.inheritance }}</td>
              <td>{{ g.variant_count }}</td>
            </tr>
          </table>
        </details>
      </template>

      <!-- Full variant catalog download -->
      <div v-if="manifest" class="catalog-download">
        <p class="note">
          The full variant catalog ({{ manifest.catalog_variant_count.toLocaleString() }}
          entries with ClinVar IDs) is embedded in the tool. Each variant links
          to its ClinVar record at
          <code>https://www.ncbi.nlm.nih.gov/clinvar/variation/{ID}/</code>
        </p>
        <button class="btn-download-small" @click="downloadCatalog">
          Download full variant catalog (CSV)
        </button>
      </div>

      <p v-if="data.liftover_note" class="note">{{ data.liftover_note }}</p>

      <p class="note audit-note">
        Every coordinate is sourced from a public database. Gene positions from
        NCBI Gene. REP coordinates from UCSC Segmental Duplications track.
        Variant catalog from ClinVar, filtered by PanelApp panel 85
        (Hereditary neuropathy). Click any source link to verify.
      </p>
    </details>
  </section>
</template>

<style scoped>
.data-sources {
  background: var(--surface);
  border-radius: var(--r-lg);
  padding: 1.25rem 1.75rem;
  margin: 1.25rem 0;
}
.data-sources .inline-heading {
  line-height: 1;
  margin: 0;
}
.inline-heading {
  display: inline;
  font-size: 1.1rem;
}
h4 {
  margin: 1rem 0 0.3rem;
  font-size: 0.95rem;
}
h5 {
  margin: 0.75rem 0 0.3rem;
  font-size: 0.85rem;
}
.detail-table {
  font-size: 0.8rem;
}
.detail-table th {
  background: var(--code-bg);
  font-weight: 600;
  padding: 0.3rem 0.5rem;
  border: 1px solid var(--border);
  text-align: left;
}
.detail-table td {
  padding: 0.3rem 0.5rem;
  border: 1px solid var(--border);
}
.audit-note {
  margin-top: 1rem;
}
.src-badge {
  font-size: 0.65rem;
  color: var(--accent);
  text-decoration: none;
  border: 1px solid var(--accent);
  border-radius: 3px;
  padding: 0 0.2rem;
  margin-left: 0.25rem;
  opacity: 0.6;
}
.src-badge:hover {
  opacity: 1;
  background: #eff6ff;
}
.src-right {
  float: right;
  margin-left: 0.5rem;
}
.catalog-download {
  margin: 1rem 0;
  padding: 0.5rem 0;
}
.btn-download-small {
  background: var(--card-bg);
  color: var(--accent);
  border: 1px solid var(--accent);
  padding: 0.3rem 0.8rem;
  border-radius: 4px;
  font-size: 0.75rem;
  cursor: pointer;
  margin-top: 0.3rem;
}
.btn-download-small:hover {
  background: #eff6ff;
}
</style>
