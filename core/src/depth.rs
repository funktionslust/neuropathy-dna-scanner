//! Global depth-ratio method.
//!
//! Computes mean read depth over the PMP22 region and the chr2 autosomal
//! control region, takes the ratio, and rounds to an integer copy-number
//! estimate. The CMT1A signal at ratio ≈ 1.5 is ~0.5 above the CN=2 baseline
//! — far outside any rounding boundary — so the simple ratio is sufficient
//! for a screening-grade result.

use crate::reader::{AlignmentSource, Region};
use crate::{AnalysisConfig, Error, Progress, Result, PROGRESS_TICK_INTERVAL};

/// Output of [`depth_ratio`]. All depth values are in reads per base;
/// `ratio` is dimensionless.
#[derive(Debug, Clone, PartialEq)]
pub struct DepthRatio {
    /// Mean read depth over the PMP22 region.
    pub pmp22_mean_depth: f64,
    /// Mean read depth over the chr2 autosomal control region.
    pub autosomal_mean_depth: f64,
    /// `pmp22_mean_depth / autosomal_mean_depth`.
    pub ratio: f64,
    /// `round(ratio * 2)`.
    pub estimated_cn: u8,
    /// Mean read depth over the RAI1 region. 0.0 if not computed.
    pub rai1_mean_depth: f64,
    /// `rai1_mean_depth / autosomal_mean_depth`. 0.0 if not computed.
    pub rai1_ratio: f64,
    /// `round(rai1_ratio * 2)`.
    pub rai1_estimated_cn: u8,
}

/// Compute the global depth ratio for a given alignment source. Walks the
/// PMP22 region and the chr2 control region via [`AlignmentSource::query_region`],
/// computes mean depth over each, and returns the ratio plus the rounded
/// copy-number estimate.
///
/// Returns [`Error::EmptyRegion`] if either region has no contributing
/// bases. This is the input-validation path for "zero mapped reads in
/// the region of interest".
pub fn depth_ratio<A: AlignmentSource>(
    source: &mut A,
    config: &AnalysisConfig,
) -> Result<DepthRatio> {
    depth_ratio_with_progress(source, config, &mut |_: Progress| {})
}

/// Progress-aware variant of [`depth_ratio`]. Emits [`Progress::WalkStart`],
/// [`Progress::WalkTick`], and [`Progress::WalkDone`] events for each
/// region walked. Used by `analyze_with_progress`.
pub fn depth_ratio_with_progress<A, F>(
    source: &mut A,
    config: &AnalysisConfig,
    on_progress: &mut F,
) -> Result<DepthRatio>
where
    A: AlignmentSource,
    F: FnMut(Progress),
{
    let pmp22_mean_depth = compute_mean_depth(
        source,
        &config.pmp22_region,
        "depth_ratio:pmp22",
        on_progress,
    )?;
    if pmp22_mean_depth == 0.0 {
        return Err(Error::EmptyRegion {
            region: config.pmp22_region.to_string(),
        });
    }

    let autosomal_mean_depth = compute_mean_depth(
        source,
        &config.chr2_control_region,
        "depth_ratio:chr2_control",
        on_progress,
    )?;
    if autosomal_mean_depth == 0.0 {
        return Err(Error::EmptyRegion {
            region: config.chr2_control_region.to_string(),
        });
    }

    let ratio = pmp22_mean_depth / autosomal_mean_depth;
    let estimated_cn = (ratio * 2.0).round().min(255.0) as u8;

    // RAI1 depth. Non-fatal if empty — RAI1 may be absent from
    // fixture BAMs or sparse byte-range slices.
    let rai1_mean_depth =
        compute_mean_depth(source, &config.rai1_region, "depth_ratio:rai1", on_progress)
            .unwrap_or(0.0);
    let rai1_ratio = if autosomal_mean_depth > 0.0 && rai1_mean_depth > 0.0 {
        rai1_mean_depth / autosomal_mean_depth
    } else {
        0.0
    };
    let rai1_estimated_cn = if rai1_ratio > 0.0 {
        (rai1_ratio * 2.0).round() as u8
    } else {
        0
    };

    Ok(DepthRatio {
        pmp22_mean_depth,
        autosomal_mean_depth,
        ratio,
        estimated_cn,
        rai1_mean_depth,
        rai1_ratio,
        rai1_estimated_cn,
    })
}

/// Compute mean depth over a region as
/// `sum(clipped_overlap_length) / region.length()`. Clipping to the region
/// boundary ensures reads that straddle the edge contribute only their
/// in-region portion to the sum. Emits progress events at `WalkStart`,
/// `WalkTick` (every [`PROGRESS_TICK_INTERVAL`] records), and `WalkDone`.
fn compute_mean_depth<A, F>(
    source: &mut A,
    region: &Region,
    phase: &'static str,
    on_progress: &mut F,
) -> Result<f64>
where
    A: AlignmentSource,
    F: FnMut(Progress),
{
    on_progress(Progress::WalkStart {
        phase,
        region: region.to_string(),
    });
    let mut total_bases: u64 = 0;
    let mut records_seen: u64 = 0;
    for record_result in source.query_region(region)? {
        let record = record_result?;
        let overlap_start = record.start.max(region.start);
        let overlap_end = record.end.min(region.end);
        if overlap_start <= overlap_end {
            total_bases += u64::from(overlap_end - overlap_start + 1);
        }
        records_seen += 1;
        if records_seen.is_multiple_of(PROGRESS_TICK_INTERVAL) {
            on_progress(Progress::WalkTick {
                phase,
                records_seen,
                last_position: record.start,
            });
        }
    }
    on_progress(Progress::WalkDone {
        phase,
        records_seen,
    });
    Ok(total_bases as f64 / f64::from(region.length()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mock::MockAlignmentSource;
    use crate::reader::AlignmentHeader;

    /// Build a mock with synthetic coverage over both regions at the given
    /// depths. `control_depth` may be 0 to skip adding chr2 records (useful
    /// for the "empty control region" error test).
    fn mock_with_coverage(
        config: &AnalysisConfig,
        pmp22_depth: u32,
        control_depth: u32,
    ) -> MockAlignmentSource {
        let mut src = MockAlignmentSource::with_coverage(&config.pmp22_region, pmp22_depth);
        if control_depth > 0 {
            src = src.add_coverage(&config.chr2_control_region, control_depth);
        }
        src
    }

    /// Empty-PMP22 variant: no pmp22 coverage, only control coverage.
    fn mock_control_only(config: &AnalysisConfig, control_depth: u32) -> MockAlignmentSource {
        MockAlignmentSource::new(AlignmentHeader::default(), Vec::new())
            .add_coverage(&config.chr2_control_region, control_depth)
    }

    #[test]
    fn depth_ratio_cn2_normal() {
        let config = AnalysisConfig::default();
        let mut mock = mock_with_coverage(&config, 30, 30);
        let result = depth_ratio(&mut mock, &config).unwrap();
        assert!(
            (result.ratio - 1.0).abs() < 0.01,
            "ratio={} (expected ~1.0)",
            result.ratio
        );
        assert_eq!(result.estimated_cn, 2);
    }

    #[test]
    fn depth_ratio_cn3_cmt1a() {
        let config = AnalysisConfig::default();
        let mut mock = mock_with_coverage(&config, 45, 30);
        let result = depth_ratio(&mut mock, &config).unwrap();
        assert!(
            (result.ratio - 1.5).abs() < 0.01,
            "ratio={} (expected ~1.5)",
            result.ratio
        );
        assert_eq!(result.estimated_cn, 3);
    }

    #[test]
    fn depth_ratio_cn1_hnpp() {
        let config = AnalysisConfig::default();
        let mut mock = mock_with_coverage(&config, 15, 30);
        let result = depth_ratio(&mut mock, &config).unwrap();
        assert!(
            (result.ratio - 0.5).abs() < 0.01,
            "ratio={} (expected ~0.5)",
            result.ratio
        );
        assert_eq!(result.estimated_cn, 1);
    }

    #[test]
    fn depth_ratio_cn4_triplication() {
        let config = AnalysisConfig::default();
        let mut mock = mock_with_coverage(&config, 60, 30);
        let result = depth_ratio(&mut mock, &config).unwrap();
        assert!(
            (result.ratio - 2.0).abs() < 0.01,
            "ratio={} (expected ~2.0)",
            result.ratio
        );
        assert_eq!(result.estimated_cn, 4);
    }

    #[test]
    fn depth_ratio_empty_pmp22_errors() {
        let config = AnalysisConfig::default();
        let mut mock = mock_control_only(&config, 30);
        let err = match depth_ratio(&mut mock, &config) {
            Ok(r) => panic!("expected EmptyRegion, got {r:?}"),
            Err(e) => e,
        };
        match err {
            Error::EmptyRegion { region } => {
                assert!(
                    region.contains("chr17"),
                    "region string should mention chr17: {region}"
                );
                assert_eq!(region, "chr17:15229779-15265326");
            }
            other => panic!("expected EmptyRegion, got {other:?}"),
        }
    }

    #[test]
    fn depth_ratio_empty_control_errors() {
        let config = AnalysisConfig::default();
        // Only pmp22 coverage, no control.
        let mut mock = mock_with_coverage(&config, 30, 0);
        let err = match depth_ratio(&mut mock, &config) {
            Ok(r) => panic!("expected EmptyRegion, got {r:?}"),
            Err(e) => e,
        };
        match err {
            Error::EmptyRegion { region } => {
                assert!(
                    region.contains("chr2"),
                    "region string should mention chr2: {region}"
                );
                assert_eq!(region, "chr2:50000000-70000000");
            }
            other => panic!("expected EmptyRegion, got {other:?}"),
        }
    }
}
