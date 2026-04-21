<script setup lang="ts">
import { computed } from "vue";
import type { RepRecord, GeneRecord, BoundaryRun } from "../types/wasm";

const props = defineProps<{
  reps: RepRecord[];
  genes: GeneRecord[];
  run: BoundaryRun;
}>();

// Visible window: extend 500 kb beyond the outermost feature
const allPositions = computed(() => {
  const positions: number[] = [];
  for (const r of props.reps) {
    positions.push(r.region.start.value, r.region.end.value);
  }
  for (const g of props.genes) {
    positions.push(g.region.start.value, g.region.end.value);
  }
  positions.push(props.run.start, props.run.end);
  return positions;
});

const windowStart = computed(
  () => Math.min(...allPositions.value) - 500_000,
);
const windowEnd = computed(
  () => Math.max(...allPositions.value) + 500_000,
);
const windowLen = computed(() => windowEnd.value - windowStart.value);

function toPercent(pos: number): number {
  return ((pos - windowStart.value) / windowLen.value) * 100;
}

function widthPercent(start: number, end: number): number {
  return ((end - start) / windowLen.value) * 100;
}

function formatMb(pos: number): string {
  return (pos / 1_000_000).toFixed(1) + " Mb";
}

// Tick marks every 1 Mb
const ticks = computed(() => {
  const result: { pos: number; label: string; pct: number }[] = [];
  const startMb = Math.ceil(windowStart.value / 1_000_000);
  const endMb = Math.floor(windowEnd.value / 1_000_000);
  for (let mb = startMb; mb <= endMb; mb++) {
    const pos = mb * 1_000_000;
    result.push({
      pos,
      label: `${mb}`,
      pct: toPercent(pos),
    });
  }
  return result;
});

const isDuplication = computed(() =>
  props.run.direction.toLowerCase().includes("dup"),
);

// Color for the detected run
const runColor = computed(() =>
  isDuplication.value ? "var(--red)" : "var(--amber)",
);
const runLabel = computed(() =>
  isDuplication.value ? "Duplication" : "Deletion",
);
</script>

<template>
  <div class="diagram">
    <div class="diagram-label">chr17</div>

    <!-- Chromosome bar -->
    <div class="chrom-bar">
      <!-- Tick marks -->
      <div
        v-for="t in ticks"
        :key="t.pos"
        class="tick"
        :style="{ left: t.pct + '%' }"
      >
        <div class="tick-line"></div>
        <div class="tick-label">{{ t.label }}</div>
      </div>

      <!-- REP elements -->
      <a
        v-for="rep in reps"
        :key="rep.name"
        :href="rep.publication_url"
        target="_blank"
        rel="noopener"
        class="feature rep-feature"
        :style="{
          left: toPercent(rep.region.start.value) + '%',
          width: Math.max(widthPercent(rep.region.start.value, rep.region.end.value), 0.8) + '%',
        }"
        :title="`${rep.name}: ${formatMb(rep.region.start.value)} - ${formatMb(rep.region.end.value)} (click for source)`"
      >
        <span class="feature-label rep-label">{{ rep.name.replace('CMT1A-REP ', '').replace('SMS-REP ', 'SMS-') }}</span>
      </a>

      <!-- Gene positions -->
      <a
        v-for="gene in genes"
        :key="gene.symbol"
        :href="gene.region.start.source_url"
        target="_blank"
        rel="noopener"
        class="feature gene-feature"
        :style="{
          left: toPercent(gene.region.start.value) + '%',
          width: Math.max(widthPercent(gene.region.start.value, gene.region.end.value), 0.6) + '%',
        }"
        :title="`${gene.symbol}: ${formatMb(gene.region.start.value)} - ${formatMb(gene.region.end.value)} (click for source)`"
      >
        <span class="feature-label gene-label">{{ gene.symbol }}</span>
      </a>

      <!-- Detected run overlay -->
      <div
        class="feature run-feature"
        :style="{
          left: toPercent(run.start) + '%',
          width: Math.max(widthPercent(run.start, run.end), 0.5) + '%',
          borderColor: runColor,
          backgroundColor: isDuplication ? 'rgba(220, 38, 38, 0.12)' : 'rgba(217, 119, 6, 0.12)',
        }"
        :title="`Your ${run.direction}: ${formatMb(run.start)} - ${formatMb(run.end)}`"
      >
        <span class="feature-label run-label" :style="{ color: runColor }">
          Your {{ runLabel.toLowerCase() }}
        </span>
      </div>
    </div>

    <!-- Legend -->
    <div class="legend">
      <span class="legend-item">
        <span class="legend-swatch rep-swatch"></span> REP elements (known breakpoint regions)
      </span>
      <span class="legend-item">
        <span class="legend-swatch gene-swatch"></span> Genes
      </span>
      <span class="legend-item">
        <span class="legend-swatch" :style="{ backgroundColor: runColor, opacity: 0.3, borderColor: runColor }"></span> Your detected region
      </span>
    </div>

    <p class="diagram-note">
      Positions in megabases (Mb) on chromosome 17. REP elements are segmental
      duplications where breakpoints typically occur. Your detected region
      should overlap with a known REP pair for a standard (recurrent)
      rearrangement.
    </p>
  </div>
</template>

<style scoped>
.diagram {
  margin: 1rem 0;
}

.diagram-label {
  font-size: 0.8rem;
  font-weight: 600;
  color: var(--text-light);
  margin-bottom: 0.25rem;
}

.chrom-bar {
  position: relative;
  height: 80px;
  background: var(--code-bg);
  border-radius: 4px;
  border: 1px solid var(--border);
  margin-bottom: 0.25rem;
}

/* Tick marks */
.tick {
  position: absolute;
  top: 0;
  height: 100%;
  pointer-events: none;
}
.tick-line {
  position: absolute;
  top: 0;
  width: 1px;
  height: 100%;
  background: var(--border);
  opacity: 0.5;
}
.tick-label {
  position: absolute;
  bottom: 2px;
  left: 2px;
  font-size: 0.6rem;
  color: var(--text-light);
  opacity: 0.6;
}

/* Features */
.feature {
  position: absolute;
  border-radius: 2px;
  min-width: 2px;
}

.rep-feature {
  top: 4px;
  height: 22px;
  background: rgba(37, 99, 235, 0.15);
  border: 1px solid var(--accent);
}

.gene-feature {
  top: 32px;
  height: 16px;
  background: rgba(22, 163, 74, 0.15);
  border: 1px solid var(--green);
}

.run-feature {
  top: 52px;
  height: 20px;
  border: 2px solid;
  border-radius: 3px;
}

.feature-label {
  position: absolute;
  top: -1px;
  left: 2px;
  font-size: 0.55rem;
  font-weight: 600;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  max-width: 100%;
  line-height: 1.1;
}

.rep-label {
  color: var(--accent);
  top: 3px;
}

.gene-label {
  color: var(--green);
  top: 1px;
  font-style: italic;
}

.run-label {
  font-size: 0.6rem;
  top: 2px;
}

/* Legend */
.legend {
  display: flex;
  gap: 1rem;
  flex-wrap: wrap;
  font-size: 0.75rem;
  color: var(--text-light);
  margin-top: 0.5rem;
}

.legend-item {
  display: flex;
  align-items: center;
  gap: 0.25rem;
}

.legend-swatch {
  display: inline-block;
  width: 14px;
  height: 10px;
  border-radius: 2px;
  border: 1px solid;
}

.rep-swatch {
  background: rgba(37, 99, 235, 0.15);
  border-color: var(--accent);
}

.gene-swatch {
  background: rgba(22, 163, 74, 0.15);
  border-color: var(--green);
}

.diagram-note {
  font-size: 0.75rem;
  color: var(--text-light);
  margin-top: 0.5rem;
  line-height: 1.4;
}
</style>
