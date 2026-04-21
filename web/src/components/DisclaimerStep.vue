<script setup lang="ts">
defineEmits<{
  accept: [];
}>();
</script>

<template>
  <section class="card disclaimer">
    <h2>Please read before continuing</h2>

    <div class="disclaimer-text">
      <p><strong>Research-grade screening tool only.</strong></p>
      <p>
        This tool is NOT a medical device, NOT a diagnostic test, and NOT a
        substitute for clinical genetic testing. It screens only for the
        specific subtypes listed. It does NOT screen for any other form of CMT,
        and it may miss atypical or mosaic cases.
      </p>
      <p>
        A positive result requires clinical confirmation by MLPA, array-CGH, or
        sequencing at a certified diagnostic laboratory before any medical
        decision.
      </p>
      <p>
        A negative result does NOT exclude CMT or other inherited neuropathies.
        Anyone with symptoms suggestive of neuropathy (foot drop, high arches,
        distal weakness, hammer toes, family history) should see a neurologist
        or clinical geneticist regardless.
      </p>
      <p>
        This tool is provided free of charge with no warranty of any kind. The
        authors assume no liability for any decisions, actions, or medical
        outcomes based on its output.
      </p>
    </div>

    <div class="action">
      <button class="btn-primary" @click="$emit('accept')">
        I understand  - this is not a medical test
      </button>
    </div>
  </section>
</template>

<style scoped>
.disclaimer-text {
  font-size: 0.9rem;
  margin: 1rem 0 1.5rem;
}
.disclaimer-text p {
  margin-bottom: 0.75rem;
}
.action {
  text-align: center;
}
</style>
