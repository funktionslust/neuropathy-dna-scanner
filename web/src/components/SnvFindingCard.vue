<script setup lang="ts">
import type { CheckedVariant, VariantCall } from "../types/wasm";

const props = defineProps<{
  variant: CheckedVariant;
}>();

function callType(call: VariantCall): string {
  if (typeof call !== "object") return "Unknown";
  if ("Heterozygous" in call) return "Heterozygous";
  if ("Homozygous" in call) return "Homozygous";
  if ("Hemizygous" in call) return "Hemizygous";
  return "Unknown";
}

function callDetail(call: VariantCall): string {
  if (typeof call !== "object") return "";
  const counts =
    "Heterozygous" in call ? call.Heterozygous :
    "Homozygous" in call ? call.Homozygous :
    null;
  if (counts) {
    const ab = (counts.allele_balance * 100).toFixed(0);
    return `${counts.ref_count} ref / ${counts.alt_count} alt reads (AB: ${ab}%)`;
  }
  if ("Hemizygous" in call) {
    const c = call.Hemizygous;
    return `${c.ref_count} ref / ${c.alt_count} alt reads (hemizygous)`;
  }
  return "";
}

function inheritanceLabel(inh: string): string {
  switch (inh) {
    case "AutosomalDominant": return "Autosomal dominant";
    case "AutosomalRecessive": return "Autosomal recessive";
    case "XLinked": return "X-linked";
    case "Both": return "AD/AR";
    default: return inh;
  }
}

function isCarrier(): boolean {
  const call = props.variant.call;
  const isHet = typeof call === "object" && "Heterozygous" in call;
  const inh = props.variant.inheritance;
  return isHet && (inh === "AutosomalRecessive" || inh === "Both");
}

function clinvarUrl(): string {
  return `https://www.ncbi.nlm.nih.gov/clinvar/variation/${props.variant.clinvar_id}/`;
}

function conditionList(): string[] {
  return (props.variant.condition ?? "")
    .split("|")
    .map((s) => s.trim())
    .filter((s) => s && s.toLowerCase() !== "not provided");
}
</script>

<template>
  <div class="snv-card" :class="{ carrier: isCarrier() }">
    <div class="snv-header">
      <span class="snv-gene">{{ variant.gene }}</span>
      <span class="snv-badge" :class="isCarrier() ? 'badge-carrier' : 'badge-found'">
        {{ isCarrier() ? "Carrier" : callType(variant.call) }}
      </span>
    </div>

    <div class="snv-variant">
      {{ variant.chrom }}:{{ variant.pos.toLocaleString() }}
      {{ variant.ref_allele }}&gt;{{ variant.alt_allele }}
      <span v-if="variant.consequence" class="snv-consequence">
        {{ variant.consequence }}
      </span>
    </div>

    <div class="snv-detail">{{ callDetail(variant.call) }}</div>

    <div v-if="conditionList().length" class="snv-condition">
      {{ conditionList().join(" / ") }}
    </div>

    <details class="snv-expand">
      <summary>Details</summary>
      <table class="detail-table">
        <tr>
          <td>ClinVar</td>
          <td>
            <a :href="clinvarUrl()" target="_blank" rel="noopener">
              #{{ variant.clinvar_id }}
            </a>
            ({{ variant.significance }}, {{ variant.review_stars }} stars)
          </td>
        </tr>
        <tr>
          <td>Inheritance</td>
          <td>{{ inheritanceLabel(variant.inheritance) }}</td>
        </tr>
        <tr v-if="isCarrier()">
          <td colspan="2" class="carrier-note">
            Carrier: this gene is autosomal recessive. A single heterozygous
            variant means you carry one copy but are unlikely to be affected
            unless a second pathogenic variant exists on the other allele.
          </td>
        </tr>
      </table>
    </details>
  </div>
</template>

<style scoped>
.snv-card {
  background: var(--red-light);
  border: 1px solid var(--red);
  border-left: 4px solid var(--red);
  border-radius: 6px;
  padding: 0.75rem 1rem;
  margin-bottom: 0.75rem;
}

.snv-card.carrier {
  background: var(--amber-light);
  border-color: var(--amber);
  border-left-color: var(--amber);
}

.snv-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  margin-bottom: 0.3rem;
}

.snv-gene {
  font-weight: 700;
  font-size: 1.05rem;
  font-style: italic;
}

.snv-badge {
  font-size: 0.75rem;
  font-weight: 600;
  padding: 0.15rem 0.5rem;
  border-radius: 4px;
}

.badge-found {
  background: var(--red);
  color: white;
}

.badge-carrier {
  background: var(--amber);
  color: white;
}

.snv-variant {
  font-family: monospace;
  font-size: 0.85rem;
  color: var(--text-light);
}

.snv-consequence {
  font-family: inherit;
  font-size: 0.75rem;
  color: var(--text-light);
  margin-left: 0.5rem;
}

.snv-detail {
  font-size: 0.8rem;
  color: var(--text-light);
  margin: 0.2rem 0;
}

.snv-condition {
  font-size: 0.85rem;
  margin-top: 0.3rem;
}

.snv-expand {
  margin-top: 0.5rem;
}

.snv-expand summary {
  font-size: 0.8rem;
  cursor: pointer;
  color: var(--text-light);
}

.carrier-note {
  font-size: 0.8rem;
  color: var(--amber);
  padding-top: 0.3rem;
}
</style>
