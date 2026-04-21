<script setup lang="ts">
import { ref, computed, watch } from "vue";
import { formatSize } from "../composables/useAnalysis";
import { useAnalysis, type DetectedFiles } from "../composables/useAnalysis";
import type { WasmModule } from "../types/wasm";

const props = defineProps<{
  wasm: WasmModule | null;
  detectedBuild: string | null;
}>();

const emit = defineEmits<{
  update: [files: DetectedFiles];
}>();

const { detectFileBuild, detectedBuild } = useAnalysis();

const genomeFile = ref<File | null>(null);
const indexFile = ref<File | null>(null);
const refFile = ref<File | null>(null);

const isCram = computed(
  () => genomeFile.value?.name.toLowerCase().endsWith(".cram") ?? false,
);
const isBam = computed(
  () => genomeFile.value?.name.toLowerCase().endsWith(".bam") ?? false,
);
const needsRef = computed(() => isCram.value);

const hasGenome = computed(() => !!genomeFile.value);
const hasIndex = computed(() => !!indexFile.value);
const hasRef = computed(() => !!refFile.value);

const indexAccept = computed(() => (isCram.value ? ".crai" : ".bai"));
const indexLabel = computed(() =>
  isCram.value ? "Select .crai file" : "Select .bai file",
);
const indexHint = computed(() =>
  isCram.value
    ? "Look for a .crai file in the same folder as your .cram"
    : "Look for a .bai file in the same folder as your .bam",
);

const buildLabel = computed(() => {
  if (!detectedBuild.value) return null;
  if (detectedBuild.value === "GRCh38") return "GRCh38 (hg38)";
  if (detectedBuild.value === "GRCh37") return "GRCh37 (hg19)";
  return detectedBuild.value;
});

const isGrch37 = computed(() => detectedBuild.value === "GRCh37");

// Emit files whenever they change
watch(
  [genomeFile, indexFile, refFile],
  () => {
    const files: DetectedFiles = {
      bam: isBam.value ? genomeFile.value : null,
      bai: isBam.value ? indexFile.value : null,
      cram: isCram.value ? genomeFile.value : null,
      crai: isCram.value ? indexFile.value : null,
      ref: refFile.value,
    };
    emit("update", files);
  },
  { deep: true },
);

async function onGenomeChange(event: Event) {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0] ?? null;
  genomeFile.value = file;
  indexFile.value = null;
  refFile.value = null;

  if (file && props.wasm) {
    const name = file.name.toLowerCase();
    await detectFileBuild(props.wasm, file, name.endsWith(".cram"));
  }
}

function onIndexChange(event: Event) {
  const input = event.target as HTMLInputElement;
  indexFile.value = input.files?.[0] ?? null;
}

function onRefChange(event: Event) {
  const input = event.target as HTMLInputElement;
  refFile.value = input.files?.[0] ?? null;
}
</script>

<template>
  <section>
    <h2>Select your files</h2>

    <!-- Step 1: Genome -->
    <div class="wizard-step" :class="{ done: hasGenome }">
      <div class="step-num">1</div>
      <div class="step-body">
        <h3>Select your genome file</h3>
        <template v-if="!hasGenome">
          <p class="note">Your BAM or CRAM file from the sequencing provider</p>
          <label class="btn-primary file-btn">
            Select .bam or .cram
            <input
              type="file"
              accept=".bam,.cram"
              hidden
              @change="onGenomeChange"
            />
          </label>
        </template>
        <span v-else class="file-ok">
          {{ genomeFile!.name }} ({{ formatSize(genomeFile!.size) }})
          <span v-if="buildLabel" class="build-tag">{{ buildLabel }}</span>
        </span>
      </div>
    </div>

    <!-- GRCh37 warning -->
    <div v-if="isGrch37 && hasGenome" class="grch37-warning">
      Your data uses <strong>GRCh37 (hg19)</strong>. This tool was primarily
      validated on GRCh38. Results should be interpreted with extra caution.
    </div>

    <!-- Step 2: Index -->
    <div v-if="hasGenome" class="wizard-step" :class="{ done: hasIndex }">
      <div class="step-num">2</div>
      <div class="step-body">
        <h3>Select the index file</h3>
        <template v-if="!hasIndex">
          <p class="note">{{ indexHint }}</p>
          <label class="btn-primary file-btn">
            {{ indexLabel }}
            <input
              type="file"
              :accept="indexAccept"
              hidden
              @change="onIndexChange"
            />
          </label>
        </template>
        <span v-else class="file-ok">{{ indexFile!.name }}</span>
      </div>
    </div>

    <!-- Step 3: Reference (CRAM only) -->
    <div
      v-if="needsRef && hasIndex"
      class="wizard-step"
      :class="{ done: hasRef }"
    >
      <div class="step-num">3</div>
      <div class="step-body">
        <h3>Select the reference genome</h3>
        <template v-if="!hasRef">
          <p class="note">
            Required for CRAM decoding. Download once, keep forever.
          </p>
          <label class="btn-primary file-btn">
            Select .fna or .fasta file
            <input
              type="file"
              accept=".fna,.fasta,.fa"
              hidden
              @change="onRefChange"
            />
          </label>
          <div class="ref-links">
            <template v-if="detectedBuild === 'GRCh38'">
              <p class="note ref-download-note">
                Your data uses <strong>GRCh38</strong>. Download the matching
                reference (~834 MB, unzip to ~3 GB):
              </p>
              <a
                href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                target="_blank"
                rel="noopener"
                class="btn-download-sm"
                >Download GRCh38 reference</a
              >
            </template>
            <template v-else-if="detectedBuild === 'GRCh37'">
              <p class="note ref-download-note">
                Your data uses <strong>GRCh37</strong>. Download the matching
                reference (~835 MB, unzip to ~3 GB):
              </p>
              <a
                href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.14_GRCh37.p13_no_alt_plus_hs38d1_analysis_set.fna.gz"
                target="_blank"
                rel="noopener"
                class="btn-download-sm"
                >Download GRCh37 reference</a
              >
            </template>
            <template v-else>
              <p class="note ref-download-note">
                Download the reference matching your data (most recent WGS uses
                GRCh38):
              </p>
              <a
                href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                target="_blank"
                rel="noopener"
                class="btn-download-sm"
                >GRCh38 / hg38</a
              >
              <a
                href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.14_GRCh37.p13_no_alt_plus_hs38d1_analysis_set.fna.gz"
                target="_blank"
                rel="noopener"
                class="btn-download-sm"
                >GRCh37 / hg19</a
              >
            </template>
          </div>
        </template>
        <span v-else class="file-ok">
          {{ refFile!.name }} ({{ formatSize(refFile!.size) }})
        </span>
      </div>
    </div>
  </section>
</template>

<style scoped>
.wizard-step {
  display: flex;
  gap: 1rem;
  margin-bottom: 1rem;
  padding: 1rem;
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-radius: 8px;
}
.wizard-step.done {
  opacity: 0.7;
}
.step-num {
  flex-shrink: 0;
  display: inline-block;
  width: 28px;
  height: 28px;
  line-height: 28px;
  text-align: center;
  background: var(--accent);
  color: #fff;
  border-radius: 50%;
  font-weight: 700;
  font-size: 0.9rem;
}
.step-body {
  flex: 1;
}
.step-body h3 {
  margin: 0 0 0.3rem;
  font-size: 1rem;
}
.file-btn {
  display: inline-block;
  cursor: pointer;
  font-size: 0.9rem;
  padding: 0.5rem 1.2rem;
  margin: 0.3rem 0;
}
.file-ok {
  color: var(--green);
  font-weight: 600;
  font-size: 0.9rem;
}
.build-tag {
  background: var(--code-bg);
  padding: 0.1rem 0.5rem;
  border-radius: 4px;
  font-size: 0.8rem;
  font-weight: normal;
  color: var(--text-light);
  margin-left: 0.5rem;
}
.ref-links {
  margin: 0.5rem 0;
}
.ref-download-note {
  margin-top: 0.5rem;
}
.btn-download-sm {
  display: inline-block;
  background: var(--card-bg);
  color: var(--accent);
  border: 1px solid var(--accent);
  padding: 0.25rem 0.6rem;
  border-radius: 4px;
  font-size: 0.75rem;
  cursor: pointer;
  text-decoration: none;
  margin-right: 0.5rem;
}
.btn-download-sm:hover {
  background: #eff6ff;
}
.grch37-warning {
  background: var(--amber-light);
  border: 1px solid var(--amber);
  border-radius: 6px;
  padding: 0.5rem 1rem;
  font-size: 0.85rem;
  margin-bottom: 1rem;
}
</style>
