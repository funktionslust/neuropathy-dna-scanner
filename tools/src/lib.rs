//! Shared utilities for neuropathy-dna-scanner development tools.
//!
//! Hosts the `synth_bam` fixture generator, the `xtask` runner, the
//! catalog builder, and the integration-fixture builder. Not a runtime
//! dependency: only `core`'s dev-dependencies pull this in.

pub mod synth;
