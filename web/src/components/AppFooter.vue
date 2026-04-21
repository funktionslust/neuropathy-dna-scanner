<script setup lang="ts">
import { computed } from "vue";
import { useWasm } from "../composables/useWasm";

defineEmits<{
  showLegal: [];
}>();

const { wasmVersion } = useWasm();
const releaseUrl = computed(() =>
  wasmVersion.value
    ? `https://github.com/funktionslust/neuropathy-dna-scanner/releases/tag/v${wasmVersion.value}`
    : "https://github.com/funktionslust/neuropathy-dna-scanner/releases",
);
</script>

<template>
  <footer class="footer">
    <p>
      Neuropathy DNA Scanner is free, open-source software (MIT licence).
      <a
        href="https://github.com/funktionslust/neuropathy-dna-scanner"
        target="_blank"
        rel="noopener"
        >Source code on GitHub</a
      >
      &middot;
      <a :href="releaseUrl" target="_blank" rel="noopener">
        Release notes<span v-if="wasmVersion"> (v{{ wasmVersion }})</span>
      </a>.
    </p>
    <p>
      Your genetic data is processed locally in your browser. No data is
      transmitted, stored, or logged.
    </p>
    <p class="credits">
      Built by Wolfgang Stark /
      <a href="https://funktionslust.com" target="_blank" rel="noopener">Funktionslust GmbH</a>
      &middot;
      <a href="mailto:info@funktionslust.digital">info@funktionslust.digital</a>
      &middot;
      <a href="#" @click.prevent="$emit('showLegal')">Legal &amp; Privacy</a>
    </p>
  </footer>
</template>

<style scoped>
.footer {
  margin-top: 3rem;
  padding-top: 1rem;
  border-top: 1px solid var(--border);
  font-size: 0.8rem;
  color: var(--text-light);
  text-align: center;
}
.footer p {
  margin: 0.2rem 0;
}
</style>
