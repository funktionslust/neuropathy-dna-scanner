# syntax=docker/dockerfile:1.7

# ---------- Stage 1: build the WebAssembly bundle ----------
# `rust:1-bookworm` tracks the latest stable 1.x release. Unpinning from
# 1.87 so that wasm-pack's toolchain requirement (>=1.89 as of 0.14) is met.
FROM rust:1-bookworm AS wasm-build

RUN rustup target add wasm32-unknown-unknown \
    && cargo install wasm-pack --locked --version 0.14.0

WORKDIR /src
COPY Cargo.toml Cargo.lock rustfmt.toml ./
COPY core/ core/
COPY tools/ tools/
COPY wasm/ wasm/
COPY vendor/ vendor/

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/src/target \
    wasm-pack build --target web --release wasm \
    && cp -r wasm/pkg /tmp/wasm-pkg

# ---------- Stage 2: build the Vue frontend ----------
FROM node:22-bookworm-slim AS web-build

WORKDIR /src/web
COPY web/package.json web/package-lock.json ./
RUN --mount=type=cache,target=/root/.npm \
    npm ci --no-audit --no-fund

COPY web/ ./
COPY --from=wasm-build /tmp/wasm-pkg ./dist/pkg
RUN npm run build

# ---------- Stage 3: serve the static bundle ----------
FROM caddy:2-alpine

COPY deploy/Caddyfile /etc/caddy/Caddyfile
COPY --from=web-build /src/web/dist /srv

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=3s --retries=3 \
    CMD wget --spider -q http://127.0.0.1:8080/ || exit 1
