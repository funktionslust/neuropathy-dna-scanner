//! Shared helpers for fixture-based integration tests.

use nds_core::reader::Region;
use nds_core::AnalysisConfig;

/// Returns an `AnalysisConfig` with fixture-sized regions matching the
/// coverage specs from `tools/src/synth.rs::fixture_specs()`.
///
/// Both this function and `fixture_specs()` use the constants from
/// `nds_tools::synth::fixture_regions` as the single source
/// of truth. No manual "keep in sync" needed.
///
/// Smaller than production defaults to keep committed BAM fixtures under
/// 20 MB total. Uses the real PMP22 coordinates but shrinks the chr2
/// control (100 kb) and boundary scan (550 kb) regions.
pub fn fixture_config() -> AnalysisConfig {
    use nds_tools::synth::fixture_regions::*;

    AnalysisConfig {
        pmp22_region: Region::new("chr17", PMP22.0, PMP22.1),
        chr2_control_region: Region::new("chr2", CHR2_CONTROL.0, CHR2_CONTROL.1),
        boundary_scan_region: Region::new("chr17", BOUNDARY.0, BOUNDARY.1),
        min_classical_run_length: MIN_CLASSICAL_RUN_LENGTH,
        max_classical_run_length: MAX_CLASSICAL_RUN_LENGTH,
        ..AnalysisConfig::default()
    }
}
