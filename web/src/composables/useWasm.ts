import { ref, shallowRef } from "vue";
import type { WasmModule, ReferenceDataset } from "../types/wasm";

const wasm = shallowRef<WasmModule | null>(null);
const loading = ref(false);
const error = ref<string | null>(null);
const wasmVersion = ref("");

let loadPromise: Promise<WasmModule> | null = null;

async function loadWasm(): Promise<WasmModule> {
  if (wasm.value) return wasm.value;
  if (loadPromise) return loadPromise;

  loading.value = true;
  error.value = null;

  loadPromise = (async () => {
    try {
      // Dynamic import of WASM glue code served from public/pkg/.
      // The path is resolved at runtime, not by the bundler.
      const mod = await import(
        /* @vite-ignore */ new URL("/pkg/nds_wasm.js", location.origin).href
      );
      await mod.default();
      const w = mod as WasmModule;
      wasm.value = w;
      wasmVersion.value = w.version();
      return w;
    } catch (e) {
      error.value = e instanceof Error ? e.message : String(e);
      throw e;
    } finally {
      loading.value = false;
    }
  })();

  return loadPromise;
}

function getReferenceData(build: string): ReferenceDataset | null {
  if (!wasm.value) return null;
  try {
    const json = wasm.value.reference_data_json(build);
    return JSON.parse(json) as ReferenceDataset;
  } catch {
    return null;
  }
}

function detectBuild(headerBytes: Uint8Array, isCram: boolean): string | null {
  if (!wasm.value) return null;
  try {
    return wasm.value.detect_build(headerBytes, isCram);
  } catch {
    return null;
  }
}

export function useWasm() {
  return {
    wasm,
    loading,
    error,
    wasmVersion,
    loadWasm,
    getReferenceData,
    detectBuild,
  };
}
