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
      <img src="/favicon.svg" alt="" class="hero-logo" width="96" height="96" />
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
    <div v-if="currentStep === 'files'" class="results-container">
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
    </div>

    <!-- Step: Results (progress + result) -->
    <template v-if="currentStep === 'result'">
      <div class="results-container">
        <AnalysisProgress :logs="logs" :state="state" :progress="progress" />
      </div>
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
/* Neuropathy DNA Scanner - quiet, medical-grade aesthetic.
   Palette uses a slate scale (neutrals) + a single indigo accent. */

:root {
  /* Neutrals - tinted periwinkle body so white panels float off it */
  --bg:          #e4e9ff;
  --surface:     #ffffff;
  --surface-2:   #f8fafc;  /* very light slate for hovers / code wells */
  --text:        #0b1120;  /* near-black, slightly warmer than pure black */
  --text-muted:  #475569;  /* slate-600 */
  --text-soft:   #64748b;  /* slate-500 */
  --border:      #e5e7eb;  /* neutral-200 */
  --border-soft: #f1f5f9;
  --code-bg:     #f8fafc;

  /* Single accent scale (indigo/blue) */
  --accent:       #2563eb;
  --accent-hover: #1d4ed8;
  --accent-soft:  #eff6ff;
  --accent-ring:  rgba(37, 99, 235, 0.18);

  /* Status colours - muted so the UI never feels alarmist */
  --green:       #059669;
  --green-soft:  #ecfdf5;
  --amber:       #b45309;
  --amber-soft:  #fffbeb;
  --red:         #b91c1c;
  --red-soft:    #fef2f2;

  /* Card + legacy aliases (for existing class names) */
  --card-bg:    var(--surface);
  --text-light: var(--text-soft);
  --green-light: var(--green-soft);
  --amber-light: var(--amber-soft);
  --red-light:   var(--red-soft);

  /* Elevation */
  --shadow-sm: 0 1px 2px rgba(15, 23, 42, 0.04);
  --shadow:    0 1px 3px rgba(15, 23, 42, 0.06), 0 1px 2px rgba(15, 23, 42, 0.04);
  --shadow-md: 0 4px 12px rgba(15, 23, 42, 0.06), 0 2px 4px rgba(15, 23, 42, 0.04);

  /* Radii */
  --r-sm: 6px;
  --r:    10px;
  --r-lg: 14px;
}

* {
  box-sizing: border-box;
  margin: 0;
  padding: 0;
}

html {
  -webkit-text-size-adjust: 100%;
  -moz-text-size-adjust: 100%;
  text-size-adjust: 100%;
}

body {
  font-family:
    "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
    "Helvetica Neue", Arial, sans-serif;
  font-feature-settings: "cv11", "ss01";
  color: var(--text);
  /* Diagonal wash: top-left stays light, bottom-right deepens into a
     richer periwinkle. Fixed attachment keeps the gradient stable as
     the page scrolls. */
  background:
    linear-gradient(135deg, #eef1ff 0%, #e4e9ff 45%, #ccd4ff 100%) fixed;
  min-height: 100vh;
  line-height: 1.55;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  font-variant-numeric: tabular-nums;
}

.app-container {
  max-width: 760px;
  margin: 0 auto;
  padding: 0 1.5rem 4rem;
}

h1, h2, h3, h4 {
  line-height: 1.25;
  letter-spacing: -0.015em;
  font-weight: 650;
}
h2 {
  margin: 2.25rem 0 0.75rem;
  font-size: 1.25rem;
}
h3 {
  font-size: 1rem;
  margin: 0.75rem 0 0.35rem;
  font-weight: 600;
}
h4 {
  font-size: 0.9rem;
  font-weight: 600;
  color: var(--text-muted);
  text-transform: none;
}

p { margin: 0; }

a {
  color: var(--accent);
  text-decoration: none;
  border-bottom: 1px solid transparent;
  transition: border-color 120ms ease;
}
a:hover { border-bottom-color: currentColor; }
a:visited { color: var(--accent); }

:focus-visible {
  outline: 2px solid var(--accent);
  outline-offset: 2px;
  border-radius: 4px;
}

/* Uniform chevron on every <summary> (details/summary expanders). */
summary {
  display: flex;
  align-items: center;
  gap: 0.6rem;
  cursor: pointer;
  list-style: none;
  user-select: none;
}
summary::-webkit-details-marker { display: none; }
summary::marker { content: ""; }
summary::before {
  content: "";
  flex-shrink: 0;
  display: inline-block;
  width: 14px;
  height: 14px;
  background-image: url('data:image/svg+xml;utf8,<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 16 16" fill="none" stroke="%2364748b" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"><path d="m6 4 4 4-4 4"/></svg>');
  background-repeat: no-repeat;
  background-position: center;
  background-size: 100% 100%;
  transition: transform 180ms ease;
}
details[open] > summary::before {
  transform: rotate(90deg);
}

.note {
  font-size: 0.825rem;
  color: var(--text-soft);
  margin-top: 0.35rem;
}

/* Disclaimer: prominent but not a bootstrap alert. Typography-led -
   the bold first line is the signal; no coloured stripe, no icon,
   no hierarchy of borders. */
.disclaimer-banner {
  background: #fef8e6;
  border: none;
  border-radius: var(--r-lg);
  padding: 1.25rem 1.5rem;
  margin: 1.5rem 0 0.5rem;
  font-size: 0.9375rem;
  line-height: 1.55;
  color: #713f12;
  font-weight: 450;
}
.disclaimer-banner strong {
  display: block;
  font-weight: 700;
  font-size: 1.05rem;
  color: #422006;
  margin-bottom: 0.4rem;
  letter-spacing: -0.015em;
}

/* Hero */
.hero {
  text-align: center;
  padding: 2.5rem 0 1rem;
}
.hero-logo {
  display: block;
  margin: 0 auto 1rem;
  filter: drop-shadow(0 6px 14px rgba(37, 99, 235, 0.18));
}
.hero h1 {
  font-size: 2rem;
  font-weight: 700;
  letter-spacing: -0.025em;
  color: var(--text);
}
.version {
  font-size: 0.75rem;
  color: var(--text-soft);
  font-weight: 500;
  font-feature-settings: "tnum";
  margin-left: 0.35rem;
  vertical-align: 3px;
  letter-spacing: 0;
}
.tagline {
  font-size: 1.0625rem;
  color: var(--text-muted);
  margin: 0.6rem 0 1rem;
  line-height: 1.5;
}
.privacy-badge {
  display: inline-flex;
  align-items: center;
  gap: 0.4rem;
  background: var(--green-soft);
  border: 1px solid #a7f3d0;
  color: #065f46;
  padding: 0.3rem 0.8rem;
  border-radius: 999px;
  font-size: 0.8rem;
  font-weight: 600;
  margin: 0.25rem 0 0.75rem;
}
.privacy-badge::before {
  content: "";
  width: 6px;
  height: 6px;
  border-radius: 50%;
  background: var(--green);
  box-shadow: 0 0 0 3px rgba(5, 150, 105, 0.18);
}
.contact-line {
  font-size: 0.825rem;
  color: var(--text-soft);
  margin: 0;
}
.contact-line a { color: var(--accent); }

/* Buttons */
.btn-primary,
.btn-secondary,
.btn-download {
  font-family: inherit;
  font-weight: 550;
  border-radius: var(--r);
  cursor: pointer;
  transition:
    background 120ms ease,
    border-color 120ms ease,
    box-shadow 120ms ease,
    transform 80ms ease;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 0.4rem;
  text-decoration: none;
  line-height: 1;
}
.btn-primary {
  background: var(--accent);
  color: #fff;
  border: 1px solid transparent;
  padding: 0.7rem 1.4rem;
  font-size: 0.95rem;
  box-shadow: var(--shadow-sm), inset 0 1px 0 rgba(255, 255, 255, 0.12);
  margin: 0.5rem 0;
}
.btn-primary:disabled {
  opacity: 0.45;
  cursor: not-allowed;
  box-shadow: none;
}
.btn-primary:hover:not(:disabled) {
  background: var(--accent-hover);
  box-shadow: var(--shadow);
}
.btn-primary:active:not(:disabled) { transform: translateY(1px); }

.btn-secondary,
.btn-download {
  background: var(--surface);
  color: var(--text);
  border: 1px solid var(--border);
  padding: 0.55rem 1.1rem;
  font-size: 0.875rem;
  box-shadow: var(--shadow-sm);
}
.btn-secondary:hover,
.btn-download:hover {
  background: var(--surface-2);
  border-color: #cbd5e1;
}

/* Surfaces - white panel floating on the tinted body, no border,
   no shadow. The body colour does the separation work. */
.card,
.results-container {
  background: var(--surface);
  border: none;
  border-radius: var(--r-lg);
  padding: 1.75rem;
  margin: 1.25rem 0;
  box-shadow: none;
}
/* Kill top margin on the first heading of a container, no matter how
   many wrappers (section, component roots) sit between it and the
   container. */
.card > *:first-child,
.results-container > *:first-child,
.card > *:first-child > *:first-child,
.results-container > *:first-child > *:first-child,
.card > *:first-child > *:first-child > *:first-child,
.results-container > *:first-child > *:first-child > *:first-child {
  margin-top: 0;
}
.card > *:last-child,
.results-container > *:last-child { margin-bottom: 0; }

.results-container > :deep(.result-card) {
  border: none;
  padding: 0;
  margin: 0 0 1rem;
  box-shadow: none;
  background: transparent;
}

/* Error surface */
.error-box {
  background: var(--red-soft);
  border: 1px solid #fecaca;
  border-left: 3px solid var(--red);
  border-radius: var(--r);
  padding: 0.75rem 1rem;
  margin: 1rem 0;
  color: #991b1b;
  font-size: 0.875rem;
}

.start-over { text-align: center; margin: 1.5rem 0 0; }

/* Tables */
.detail-table {
  width: 100%;
  border-collapse: separate;
  border-spacing: 0;
  margin: 0.5rem 0;
  font-size: 0.875rem;
}
.detail-table td,
.detail-table th {
  padding: 0.5rem 0.65rem;
  border-bottom: 1px solid var(--border-soft);
  text-align: left;
  vertical-align: top;
}
.detail-table tr:last-child td,
.detail-table tr:last-child th {
  border-bottom: none;
}
.detail-table td:first-child {
  color: var(--text-soft);
  width: 40%;
  font-size: 0.8125rem;
}
.detail-table th {
  background: var(--surface-2);
  font-weight: 600;
  border: 1px solid var(--border);
  color: var(--text-muted);
  font-size: 0.8125rem;
}

code {
  font-family: "SF Mono", ui-monospace, SFMono-Regular, Menlo, Monaco,
    "Cascadia Mono", Consolas, monospace;
  font-size: 0.85em;
  background: var(--code-bg);
  padding: 0.05rem 0.35rem;
  border-radius: 4px;
  color: var(--text-muted);
}

/* Legal overlay */
.legal-overlay {
  position: fixed;
  inset: 0;
  background: rgba(15, 23, 42, 0.55);
  backdrop-filter: blur(6px);
  z-index: 1000;
  overflow-y: auto;
  padding: 2rem 1rem;
  animation: overlay-in 180ms ease;
}
@keyframes overlay-in {
  from { opacity: 0; backdrop-filter: blur(0); }
  to   { opacity: 1; backdrop-filter: blur(6px); }
}

/* Responsive */
@media (max-width: 600px) {
  .app-container { padding: 0 1rem 2.5rem; }
  .hero          { padding: 1.5rem 0 0.5rem; }
  .hero h1       { font-size: 1.5rem; }
  .tagline       { font-size: 1rem; }
  .card,
  .results-container { padding: 1.1rem; border-radius: var(--r); }
}

@media (prefers-color-scheme: dark) {
  :root {
    --bg:          #0b1120;
    --surface:     #0f172a;
    --surface-2:   #111827;
    --text:        #e2e8f0;
    --text-muted:  #94a3b8;
    --text-soft:   #64748b;
    --border:      #1e293b;
    --border-soft: #172036;
    --code-bg:     #111827;
    --accent-soft: rgba(37, 99, 235, 0.14);
    --green-soft:  rgba(5, 150, 105, 0.12);
    --amber-soft:  rgba(180, 83, 9, 0.14);
    --red-soft:    rgba(185, 28, 28, 0.14);
    --shadow-sm:   0 1px 2px rgba(0, 0, 0, 0.35);
    --shadow:      0 1px 3px rgba(0, 0, 0, 0.5);
    --shadow-md:   0 4px 12px rgba(0, 0, 0, 0.5);
  }
  .disclaimer-banner { color: #fcd34d; }
  .disclaimer-banner strong { color: #fef3c7; }
  .privacy-badge { color: #6ee7b7; border-color: #064e3b; }
  .error-box { color: #fca5a5; border-color: #7f1d1d; }
}
</style>
