<script setup lang="ts">
import { computed, ref, watch, nextTick } from "vue";
import type { AnalysisState, ProgressInfo } from "../composables/useAnalysis";

const props = defineProps<{
  logs: string[];
  state: AnalysisState;
  progress: ProgressInfo;
}>();

const logEl = ref<HTMLPreElement | null>(null);
const logsExpanded = ref(false);

const isActive = computed(
  () => props.state === "loading" || props.state === "analyzing",
);

const isDone = computed(() => props.state === "done");
const isError = computed(() => props.state === "error");

const barWidth = computed(() => {
  if (props.progress.percent < 0) return "100%"; // indeterminate
  return `${Math.min(100, props.progress.percent)}%`;
});

const barClass = computed(() => {
  if (isError.value) return "bar-error";
  if (isDone.value) return "bar-done";
  if (props.progress.percent < 0) return "bar-pulse";
  return "";
});

watch(
  () => props.logs.length,
  async () => {
    await nextTick();
    if (logEl.value) {
      logEl.value.scrollTop = logEl.value.scrollHeight;
    }
  },
);
</script>

<template>
  <section class="progress-section">
    <!-- Phase label -->
    <div class="phase-row">
      <span v-if="isDone && progress.totalSteps > 0" class="step-tag">
        Step {{ progress.totalSteps }}/{{ progress.totalSteps }}
      </span>
      <span v-else-if="progress.step > 0 && progress.totalSteps > 0" class="step-tag">
        Step {{ progress.step }}/{{ progress.totalSteps }}
      </span>
      <span class="phase-label">{{ progress.phase }}</span>
      <span v-if="isActive && progress.percent >= 0" class="phase-pct">
        {{ Math.round(progress.percent) }}%
      </span>
      <span v-if="isDone" class="phase-done">Done</span>
      <span v-if="isError" class="phase-error">Error</span>
    </div>
    <div v-if="isActive && progress.detail" class="phase-detail">
      {{ progress.detail }}
    </div>

    <!-- Progress bar -->
    <div class="bar-track">
      <div
        class="bar-fill"
        :class="barClass"
        :style="{ width: barWidth }"
      ></div>
    </div>

    <!-- Collapsible log -->
    <details class="log-details" :open="logsExpanded || isError">
      <summary class="log-toggle" @click.prevent="logsExpanded = !logsExpanded">
        {{ logsExpanded ? 'Hide' : 'Show' }} technical log
        <span class="log-count">({{ logs.length }} lines)</span>
      </summary>
      <pre ref="logEl" class="progress-log">{{ logs.join('\n') }}</pre>
    </details>
  </section>
</template>

<style scoped>
.progress-section {
  margin: 1.5rem 0;
}

.phase-row {
  display: flex;
  align-items: baseline;
  gap: 0.5rem;
  margin-bottom: 0.5rem;
}

.phase-label {
  font-weight: 600;
  font-size: 1rem;
}

.step-tag {
  font-size: 0.75rem;
  color: var(--accent);
  background: #eff6ff;
  border: 1px solid var(--accent);
  border-radius: 3px;
  padding: 0.1rem 0.4rem;
  font-weight: 600;
}

.phase-pct {
  font-size: 0.85rem;
  color: var(--text-light);
  font-variant-numeric: tabular-nums;
}

.phase-detail {
  font-size: 0.8rem;
  color: var(--text-light);
  margin-bottom: 0.4rem;
}

.phase-done {
  color: var(--green);
  font-weight: 600;
  font-size: 0.9rem;
}

.phase-error {
  color: var(--red);
  font-weight: 600;
  font-size: 0.9rem;
}

/* Progress bar */
.bar-track {
  height: 8px;
  background: var(--code-bg);
  border-radius: 4px;
  overflow: hidden;
  margin-bottom: 0.75rem;
}

.bar-fill {
  height: 100%;
  background: var(--accent);
  border-radius: 4px;
  transition: width 0.3s ease;
}

.bar-done {
  background: var(--green);
}

.bar-error {
  background: var(--red);
}

.bar-pulse {
  animation: pulse 1.5s ease-in-out infinite;
}

@keyframes pulse {
  0%,
  100% {
    opacity: 0.6;
  }
  50% {
    opacity: 1;
  }
}

/* Log toggle */
.log-details {
  margin-top: 0.5rem;
}

.log-toggle {
  cursor: pointer;
  font-size: 0.8rem;
  color: var(--text-light);
  user-select: none;
}

.log-count {
  color: var(--text-light);
  font-size: 0.75rem;
}

.progress-log {
  background: var(--code-bg);
  padding: 0.75rem;
  border-radius: 4px;
  font-size: 0.75rem;
  max-height: 250px;
  overflow-y: auto;
  white-space: pre-wrap;
  margin-top: 0.5rem;
  color: var(--text-light);
}
</style>
