<script setup lang="ts">
defineEmits<{
  continue: [];
}>();
</script>

<template>
  <section class="card welcome">
    <h2>What is this?</h2>
    <p>
      <strong>Neuropathy DNA Scanner</strong> screens your whole-genome sequencing (WGS)
      data for inherited neuropathies. It checks for large
      duplications and deletions on chromosome 17p as well as known pathogenic
      point mutations across 165+ neuropathy-associated genes from the ClinVar database.
    </p>

    <div class="conditions-list">
      <h3>What it screens for</h3>
      <ul>
        <li>
          <strong>Structural variants</strong> - large duplications and deletions
          (CMT1A, HNPP, Potocki-Lupski, Smith-Magenis, YUHAL)
        </li>
        <li>
          <strong>Point mutations</strong> - known pathogenic SNVs across
          165+ genes including PMP22, GJB1, MPZ, MFN2, GDAP1, SH3TC2, NEFL,
          EGR2 and many more (sourced from ClinVar + Genomics England PanelApp)
        </li>
      </ul>
    </div>

    <div class="how-it-works">
      <h3>How it works</h3>
      <div class="steps">
        <div class="step-card">
          <div class="step-num">1</div>
          <h4>Get your data</h4>
          <p>
            Order a WGS test from a privacy-conscious provider. You need at
            least 15x coverage. Download your BAM or CRAM file.
          </p>
        </div>
        <div class="step-card">
          <div class="step-num">2</div>
          <h4>Select your files</h4>
          <p>
            Pick your genome file and its index. The analysis runs entirely in
            your browser  -
            <strong>nothing is uploaded anywhere</strong>.
          </p>
        </div>
        <div class="step-card">
          <div class="step-num">3</div>
          <h4>Read your result</h4>
          <p>
            The tool reports any structural variants or point mutations found.
            If positive, take the report to a clinical geneticist.
          </p>
        </div>
      </div>
    </div>

    <p class="privacy-note">
      Your genetic data is processed locally in your browser using WebAssembly.
      No data is transmitted, stored, or logged. A 45 GB genome file would take
      hours to upload anyway  - that's why we analyse it locally.
    </p>

    <div class="action">
      <button class="btn-primary" @click="$emit('continue')">
        I understand, continue
      </button>
    </div>
  </section>
</template>

<style scoped>
.welcome p {
  margin: 0.5rem 0;
}
.conditions-list {
  margin: 1rem 0;
}
.conditions-list ul {
  margin: 0.5rem 0 0 1.2rem;
  font-size: 0.9rem;
}
.conditions-list li {
  margin-bottom: 0.3rem;
}
.steps {
  display: flex;
  gap: 1rem;
  margin: 1rem 0;
  flex-wrap: wrap;
}
.step-card {
  flex: 1 1 200px;
  background: var(--code-bg);
  border-radius: 6px;
  padding: 1rem;
}
.step-num {
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
  margin-bottom: 0.5rem;
}
.step-card h4 {
  margin: 0.3rem 0;
  font-size: 0.95rem;
}
.step-card p {
  font-size: 0.85rem;
  color: var(--text-light);
}
.privacy-note {
  font-size: 0.85rem;
  color: var(--text-light);
  margin: 1rem 0;
}
.action {
  text-align: center;
  margin-top: 1.5rem;
}
</style>
