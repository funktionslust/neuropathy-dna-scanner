import { defineConfig } from "vite";
import vue from "@vitejs/plugin-vue";

export default defineConfig({
  plugins: [vue()],
  root: ".",
  publicDir: "public",
  build: {
    outDir: "dist",
    // `dist/pkg/` is populated separately (freshly-built WASM pkg dropped
    // in by the Docker build / local `wasm-pack build` + cp). Don't let
    // vite wipe it when it rebuilds the Vue bundle.
    emptyOutDir: false,
  },
  server: {
    port: 8090,
  },
});
