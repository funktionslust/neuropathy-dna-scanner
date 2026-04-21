<script setup lang="ts">
import { ref, onMounted } from "vue";
import { useWasm } from "./composables/useWasm";
import { useAnalysis, type DetectedFiles } from "./composables/useAnalysis";
import WelcomeStep from "./components/WelcomeStep.vue";
import DisclaimerStep from "./components/DisclaimerStep.vue";
import FileWizard from "./components/FileWizard.vue";
import SyndromeSelector from "./components/SyndromeSelector.vue";
import AnalysisProgress from "./components/AnalysisProgress.vue";
import UnifiedResults from "./components/UnifiedResults.vue";
import DataSources from "./components/DataSources.vue";
import AppFooter from "./components/AppFooter.vue";
import LegalPage from "./components/LegalPage.vue";

const { wasm, wasmVersion, loadWasm, error: wasmError } = useWasm();
const {
  state,
  logs,
  result,
  pileupResult,
  screeningProgress,
  progress,
  detectedBuild,
  detectedSex,
  isRunning,
  analyze,
  reset,
} = useAnalysis();

// Wizard steps: welcome -> disclaimer -> files -> result
type Step = "welcome" | "disclaimer" | "files" | "result";
const currentStep = ref<Step>("welcome");
const files = ref<DetectedFiles>({
  bam: null,
  bai: null,
  cram: null,
  crai: null,
  ref: null,
});
const disclaimerAccepted = ref(false);
const showLegal = ref(false);

onMounted(async () => {
  try {
    await loadWasm();
  } catch {
    // error is in wasmError
  }
  if (localStorage.getItem("cmt-disclaimer-accepted") === "yes") {
    disclaimerAccepted.value = true;
    currentStep.value = "files";
  }
});

function onWelcomeContinue() {
  if (disclaimerAccepted.value) {
    currentStep.value = "files";
  } else {
    currentStep.value = "disclaimer";
  }
}

function onDisclaimerAccept() {
  localStorage.setItem("cmt-disclaimer-accepted", "yes");
  disclaimerAccepted.value = true;
  currentStep.value = "files";
}

function onFilesReady(detected: DetectedFiles) {
  files.value = detected;
}

async function onAnalyze() {
  if (!wasm.value) return;
  currentStep.value = "result";
  await analyze(wasm.value, files.value);
}

function onStartOver() {
  reset();
  files.value = { bam: null, bai: null, cram: null, crai: null, ref: null };
  currentStep.value = "files";
}
</script>

<template>
  <div class="app-container">
    <!-- Persistent disclaimer banner -->
    <div class="disclaimer-banner">
      <strong>This is NOT a medical device and NOT a diagnostic test.</strong>
      It is a free research tool that gives you a first impression. It cannot
      replace a clinical geneticist. No liability is assumed for any decision
      made based on this tool's output.
    </div>

    <!-- Header -->
    <header class="hero">
      <h1>
        Neuropathy DNA Scanner
        <span class="version">{{
          wasmVersion ? `v${wasmVersion}` : "loading..."
        }}</span>
      </h1>
      <p class="tagline">
        Screen your genome for inherited neuropathies. Locally in your
        browser.
      </p>
      <div class="privacy-badge">Your genetic data never leaves your computer.</div>
      <p class="contact-line">
        Stuck, confused by a result, or found a bug?
        <a href="mailto:info@funktionslust.digital">info@funktionslust.digital</a>
      </p>
    </header>

    <div v-if="wasmError" class="error-box">
      WASM module failed to load: {{ wasmError }}
    </div>

    <!-- Step: Welcome -->
    <WelcomeStep
      v-if="currentStep === 'welcome'"
      @continue="onWelcomeContinue"
    />

    <!-- Step: Disclaimer -->
    <DisclaimerStep
      v-if="currentStep === 'disclaimer'"
      @accept="onDisclaimerAccept"
    />

    <!-- Step: File selection + syndrome selection + analyze -->
    <template v-if="currentStep === 'files'">
      <FileWizard
        :wasm="wasm"
        :detected-build="detectedBuild"
        @update="onFilesReady"
      />
      <SyndromeSelector
        :files="files"
        :is-running="isRunning"
        :wasm="wasm"
        @analyze="onAnalyze"
      />
    </template>

    <!-- Step: Results (progress + result) -->
    <template v-if="currentStep === 'result'">
      <AnalysisProgress :logs="logs" :state="state" :progress="progress" />
      <div v-if="result" class="results-container">
        <UnifiedResults
          :sv-result="result"
          :pileup-result="pileupResult"
          :screening="screeningProgress"
          :detected-build="detectedBuild"
          :detected-sex="detectedSex"
          :wasm="wasm"
        />
      </div>
      <div v-if="state === 'done' || state === 'error'" class="start-over">
        <button class="btn-secondary" @click="onStartOver">
          Analyze another file
        </button>
      </div>
    </template>

    <!-- Data sources (always visible, collapsed) -->
    <DataSources :wasm="wasm" :detected-build="detectedBuild" />

    <AppFooter @show-legal="showLegal = true" />

    <!-- Legal overlay -->
    <div v-if="showLegal" class="legal-overlay" @click.self="showLegal = false">
      <LegalPage @close="showLegal = false" />
    </div>
  </div>
</template>

<style>
/* Neuropathy DNA Scanner web: medical-grade, calm, trustworthy */

:root {
  --bg: #fafbfc;
  --text: #1a1a2e;
  --text-light: #555;
  --border: #dde;
  --accent: #2563eb;
  --green: #16a34a;
  --amber: #d97706;
  --red: #dc2626;
  --red-light: #fef2f2;
  --amber-light: #fffbeb;
  --green-light: #f0fdf4;
  --card-bg: #fff;
  --code-bg: #f3f4f6;
}

* {
  box-sizing: border-box;
  margin: 0;
  padding: 0;
}

body {
  font-family:
    -apple-system,
    BlinkMacSystemFont,
    "Segoe UI",
    Roboto,
    "Helvetica Neue",
    Arial,
    sans-serif;
  color: var(--text);
  background: var(--bg);
  line-height: 1.6;
}

.app-container {
  max-width: 800px;
  margin: 0 auto;
  padding: 0 1.5rem 3rem;
}

h1,
h2,
h3 {
  line-height: 1.3;
}
h2 {
  margin: 2rem 0 0.75rem;
  font-size: 1.3rem;
  color: var(--text);
}
h3 {
  font-size: 1.05rem;
  margin: 0.5rem 0;
}

a {
  color: var(--accent);
}
a:visited {
  color: #6d28d9;
}

.note {
  font-size: 0.85rem;
  color: var(--text-light);
  margin-top: 0.3rem;
}

/* Disclaimer banner */
.disclaimer-banner {
  background: var(--amber-light);
  border: 2px solid var(--amber);
  border-radius: 6px;
  padding: 0.75rem 1rem;
  margin: 1rem 0;
  font-size: 0.9rem;
  line-height: 1.5;
}

/* Hero */
.hero {
  text-align: center;
  padding: 1.5rem 0 0.5rem;
}
.hero h1 {
  font-size: 1.8rem;
}
.version {
  font-size: 0.9rem;
  color: var(--text-light);
  font-weight: normal;
}
.tagline {
  font-size: 1.1rem;
  margin: 0.5rem 0;
}
.privacy-badge {
  display: inline-block;
  background: var(--green-light);
  border: 1px solid var(--green);
  color: var(--green);
  padding: 0.3rem 0.8rem;
  border-radius: 20px;
  font-size: 0.85rem;
  font-weight: 600;
  margin: 0.5rem 0;
}
.contact-line {
  font-size: 0.85rem;
  color: var(--text-light);
  margin: 0.25rem 0 0;
}
.contact-line a {
  color: var(--accent);
}

/* Buttons */
.btn-primary {
  background: var(--accent);
  color: #fff;
  border: none;
  padding: 0.6rem 1.5rem;
  border-radius: 6px;
  font-size: 1rem;
  cursor: pointer;
  margin: 0.5rem 0;
}
.btn-primary:disabled {
  opacity: 0.4;
  cursor: not-allowed;
}
.btn-primary:hover:not(:disabled) {
  background: #1d4ed8;
}

.btn-secondary {
  background: var(--card-bg);
  color: var(--accent);
  border: 1px solid var(--accent);
  padding: 0.5rem 1rem;
  border-radius: 6px;
  font-size: 0.9rem;
  cursor: pointer;
}
.btn-secondary:hover {
  background: #eff6ff;
}

.btn-download {
  background: var(--card-bg);
  color: var(--accent);
  border: 1px solid var(--accent);
  padding: 0.5rem 1rem;
  border-radius: 6px;
  font-size: 0.9rem;
  cursor: pointer;
  text-decoration: none;
  display: inline-block;
}
.btn-download:hover {
  background: #eff6ff;
}

/* Cards */
.card {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 1.5rem;
  margin: 1rem 0;
}

/* Error */
.error-box {
  background: var(--red-light);
  border: 1px solid var(--red);
  border-radius: 6px;
  padding: 0.75rem 1rem;
  margin: 1rem 0;
  color: var(--red);
}

.results-container {
  background: var(--card-bg);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 1.5rem;
  margin: 1.5rem 0;
}

.results-container > :deep(.result-card) {
  border: none;
  padding: 0;
  margin: 0 0 1rem;
}

.start-over {
  text-align: center;
  margin: 1.5rem 0;
}

/* Detail table */
.detail-table {
  width: 100%;
  border-collapse: collapse;
  margin: 0.5rem 0;
  font-size: 0.9rem;
}
.detail-table td,
.detail-table th {
  padding: 0.3rem 0.5rem;
  border-bottom: 1px solid var(--border);
  text-align: left;
}
.detail-table td:first-child {
  color: var(--text-light);
  width: 40%;
}
.detail-table th {
  background: var(--code-bg);
  font-weight: 600;
  border: 1px solid var(--border);
}

/* Legal overlay */
.legal-overlay {
  position: fixed;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: rgba(0, 0, 0, 0.5);
  z-index: 1000;
  overflow-y: auto;
  padding: 2rem 1rem;
}

/* Responsive */
@media (max-width: 600px) {
  .app-container {
    padding: 0 1rem 2rem;
  }
  .hero h1 {
    font-size: 1.4rem;
  }
}
</style>
