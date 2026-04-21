<script setup lang="ts">
import { computed } from "vue";
import type { DetectedFiles } from "../composables/useAnalysis";
import type { WasmModule } from "../types/wasm";

const props = defineProps<{
  files: DetectedFiles;
  isRunning: boolean;
  wasm: WasmModule | null;
}>();

const emit = defineEmits<{
  analyze: [];
}>();

const canAnalyze = computed(() => {
  const f = props.files;
  const hasBam = !!(f.bam && f.bai);
  const hasCram = !!(f.cram && f.crai && f.ref);
  return (hasBam || hasCram) && !!props.wasm;
});

const isCram = computed(() => !!(props.files.cram && props.files.crai));

const whyDisabled = computed(() => {
  const f = props.files;
  if (!f.bam && !f.cram) return "Select a BAM or CRAM file first.";
  if (f.bam && !f.bai) return "BAM found but the .bai index is missing.";
  if (f.cram && !f.crai) return "CRAM found but the .crai index is missing.";
  if (isCram.value && !f.ref)
    return "CRAM needs a reference FASTA. Select it above.";
  if (!props.wasm) return "WASM module is still loading...";
  return "";
});
</script>

<template>
  <div v-if="canAnalyze" class="analyse-section">
    <div class="wizard-step">
      <div class="step-num">{{ isCram ? '4' : '3' }}</div>
      <div class="step-body">
        <h3>Analyse</h3>
        <p class="note">
          All supported syndromes will be screened. Results appear in seconds.
        </p>
        <button
          class="btn-primary"
          :disabled="!canAnalyze || isRunning"
          @click="emit('analyze')"
        >
          {{ isRunning ? "Analysing..." : "Analyse" }}
        </button>
      </div>
    </div>
  </div>
</template>

<style scoped>
.analyse-section {
  margin-top: 0.5rem;
}
.wizard-step {
  display: flex;
  gap: 1rem;
  padding: 0.75rem 0;
  border-top: 1px solid var(--border-soft);
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
</style>
