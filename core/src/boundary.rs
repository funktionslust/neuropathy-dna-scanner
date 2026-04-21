//! Sliding-window boundary scan.
//!
//! Walks the boundary-scan region (default `chr17:13,500,000-21,000,000`) in
//! 10 kb windows, normalises each window's depth against the autosomal
//! baseline, and reports the longest contiguous run of non-normal
//! (duplicated or deleted) windows.
//!
//! The region extends to 21 Mb to cover the SMS-REP proximal element at
//! ~20.3 Mb, enabling detection of PTLS, SMS, and YUHAL syndromes in
//! addition to CMT1A/HNPP.
//!
//! This is the location-and-extent complement to the global `depth_ratio`
//! method. `analyze` runs both against the same source; `interpret` combines
//! both verdicts into a final classification.

use crate::reader::{AlignmentSource, Region};
use crate::{AnalysisConfig, Error, Progress, Result, PROGRESS_TICK_INTERVAL};

/// Window size in bases.
const WINDOW_SIZE: u32 = 10_000;

/// Direction of a boundary run - duplication (elevated depth) or deletion
/// (depressed depth).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryDirection {
    Duplication,
    Deletion,
}

/// Per-window depth measurements.
#[derive(Debug, Clone, PartialEq)]
pub struct WindowDepth {
    /// 1-based inclusive window start.
    pub start: u32,
    /// 1-based inclusive window end.
    pub end: u32,
    /// Mean read depth across the window, in reads per base.
    pub raw_depth: f64,
    /// `raw_depth / autosomal_baseline`, dimensionless.
    pub normalized_depth: f64,
    /// Fraction of reads with MAPQ >= 20 (unique mappers). Regions with
    /// low unique-mapper fraction contain paralog-rich sequence where
    /// short-read depth is unreliable for CNV calling.
    pub unique_mapper_fraction: f64,
}

impl WindowDepth {
    fn direction(
        &self,
        dup_threshold: f64,
        del_threshold: f64,
        min_mapq_fraction: f64,
    ) -> Option<BoundaryDirection> {
        // Low-mappability windows are treated as "no evidence" regardless
        // of apparent depth. Paralog-rich regions produce unreliable depth
        // signals that cause false-positive CNV calls.
        if self.unique_mapper_fraction < min_mapq_fraction {
            return None;
        }
        if self.normalized_depth > dup_threshold {
            Some(BoundaryDirection::Duplication)
        } else if self.normalized_depth < del_threshold {
            Some(BoundaryDirection::Deletion)
        } else {
            None
        }
    }
}

/// A contiguous run of same-direction non-normal windows.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryRun {
    /// 1-based inclusive start coordinate of the run.
    pub start: u32,
    /// 1-based inclusive end coordinate of the run.
    pub end: u32,
    /// `end - start + 1`, in bases.
    pub length: u32,
    /// Mean normalised depth across all windows in the run.
    pub mean_normalized_depth: f64,
    /// Whether the run is a duplication or a deletion.
    pub direction: BoundaryDirection,
    /// Number of below-threshold windows bridged inside this run.
    /// Zero means a clean, uninterrupted CNV signal. Non-zero
    /// indicates segdup-edge dips that were bridged by the gap
    /// tolerance setting.
    pub bridged_gaps: u32,
}

/// Output of [`boundary_scan`]. Contains every window's depth plus all
/// detected runs of non-normal windows, sorted by length descending.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryScanResult {
    pub windows: Vec<WindowDepth>,
    /// All runs above `min_reportable_run_length`, sorted by length desc.
    /// The first element (if any) is the longest run.
    pub runs: Vec<BoundaryRun>,
}

impl BoundaryScanResult {
    /// The longest run, if any. Convenience accessor for code that only
    /// needs the primary finding.
    pub fn longest_run(&self) -> Option<&BoundaryRun> {
        self.runs.first()
    }
}

/// Scan the boundary-scan region in 10 kb windows and detect the longest
/// contiguous run of non-normal windows.
///
/// `autosomal_baseline` must be positive (typically the
/// `autosomal_mean_depth` from a successful [`crate::depth_ratio`] call). A
/// non-positive baseline returns [`Error::Internal`] rather than propagating
/// a NaN.
pub fn boundary_scan<A: AlignmentSource>(
    source: &mut A,
    autosomal_baseline: f64,
    config: &AnalysisConfig,
) -> Result<BoundaryScanResult> {
    boundary_scan_with_progress(source, autosomal_baseline, config, &mut |_: Progress| {})
}

/// Progress-aware variant of [`boundary_scan`]. Emits [`Progress::WalkStart`],
/// [`Progress::WalkTick`], and [`Progress::WalkDone`] for the single chr17
/// region walk. Used by `analyze_with_progress`.
pub fn boundary_scan_with_progress<A, F>(
    source: &mut A,
    autosomal_baseline: f64,
    config: &AnalysisConfig,
    on_progress: &mut F,
) -> Result<BoundaryScanResult>
where
    A: AlignmentSource,
    F: FnMut(Progress),
{
    if !autosomal_baseline.is_finite() || autosomal_baseline <= 0.0 {
        return Err(Error::Internal(
            "autosomal baseline must be positive".into(),
        ));
    }
    let windows = compute_window_depths(
        source,
        &config.boundary_scan_region,
        autosomal_baseline,
        on_progress,
    )?;
    let mut runs = find_all_runs(
        &windows,
        config.max_run_gap_windows,
        config.duplication_threshold,
        config.deletion_threshold,
        config.min_unique_mapper_fraction,
    );
    runs.retain(|run| run.length >= config.min_reportable_run_length);
    Ok(BoundaryScanResult { windows, runs })
}

/// Walk the boundary-scan region in a single `query_region` call, attribute
/// each record's clipped overlap to the right window(s), and return the
/// resulting window-depth vector. Emits progress events during the walk.
fn compute_window_depths<A, F>(
    source: &mut A,
    region: &Region,
    autosomal_baseline: f64,
    on_progress: &mut F,
) -> Result<Vec<WindowDepth>>
where
    A: AlignmentSource,
    F: FnMut(Progress),
{
    let num_windows = region.length() / WINDOW_SIZE;
    if num_windows == 0 {
        return Ok(Vec::new());
    }
    // Floor division drops the trailing partial window. `effective_end` is
    // the last base covered by full windows.
    let effective_end = region.start + num_windows * WINDOW_SIZE - 1;
    let n = num_windows as usize;

    let mut total_bases = vec![0u64; n];
    let mut total_reads = vec![0u64; n]; // all reads per window
    let mut unique_reads = vec![0u64; n]; // reads with MAPQ >= 20

    const PHASE: &str = "boundary_scan:chr17";
    const MAPQ_UNIQUE: u8 = 20;
    on_progress(Progress::WalkStart {
        phase: PHASE,
        region: region.to_string(),
    });
    let mut records_seen: u64 = 0;
    for record_result in source.query_region(region)? {
        let record = record_result?;
        let rstart = record.start.max(region.start);
        let rend = record.end.min(effective_end);
        if rstart > rend {
            records_seen += 1;
            continue;
        }
        let is_unique = record.mapq >= MAPQ_UNIQUE;
        let first_window = (rstart - region.start) / WINDOW_SIZE;
        let last_window = (rend - region.start) / WINDOW_SIZE;
        for w in first_window..=last_window {
            let wi = w as usize;
            let window_start = region.start + w * WINDOW_SIZE;
            let window_end = window_start + WINDOW_SIZE - 1;
            let overlap_start = rstart.max(window_start);
            let overlap_end = rend.min(window_end);
            if overlap_start <= overlap_end {
                total_bases[wi] += u64::from(overlap_end - overlap_start + 1);
                total_reads[wi] += 1;
                if is_unique {
                    unique_reads[wi] += 1;
                }
            }
        }
        records_seen += 1;
        if records_seen.is_multiple_of(PROGRESS_TICK_INTERVAL) {
            on_progress(Progress::WalkTick {
                phase: PHASE,
                records_seen,
                last_position: record.start,
            });
        }
    }
    on_progress(Progress::WalkDone {
        phase: PHASE,
        records_seen,
    });

    let mut windows = Vec::with_capacity(n);
    for i in 0..n {
        let start = region.start + (i as u32) * WINDOW_SIZE;
        let end = start + WINDOW_SIZE - 1;
        let raw_depth = total_bases[i] as f64 / f64::from(WINDOW_SIZE);
        let normalized_depth = raw_depth / autosomal_baseline;
        let unique_mapper_fraction = if total_reads[i] > 0 {
            unique_reads[i] as f64 / total_reads[i] as f64
        } else {
            1.0 // no reads = no evidence either way, don't penalise
        };
        windows.push(WindowDepth {
            start,
            end,
            raw_depth,
            normalized_depth,
            unique_mapper_fraction,
        });
    }
    Ok(windows)
}

/// Walk the window list and return ALL runs of same-direction non-normal
/// windows, sorted by length descending. Allows up to `max_gap_windows`
/// consecutive below-threshold windows inside a run (to bridge CMT1A-REP
/// segdup dips). A run's boundaries are the outermost non-normal windows,
/// i.e. trailing gaps are trimmed and never contribute to the length or
/// mean. Gap windows inside the run DO contribute to the mean, which
/// honestly reflects the average normalised depth across the CNV event.
fn find_all_runs(
    windows: &[WindowDepth],
    max_gap_windows: usize,
    dup_threshold: f64,
    del_threshold: f64,
    min_mapq_fraction: f64,
) -> Vec<BoundaryRun> {
    #[derive(Clone, Copy)]
    struct RunState {
        first_idx: usize,
        last_hit_idx: usize,
        direction: BoundaryDirection,
        gaps: u32,
    }

    let mut completed: Vec<RunState> = Vec::new();
    let mut current: Option<RunState> = None;

    for (i, w) in windows.iter().enumerate() {
        let dir = w.direction(dup_threshold, del_threshold, min_mapq_fraction);
        match (current, dir) {
            (None, None) => {}
            (None, Some(d)) => {
                current = Some(RunState {
                    first_idx: i,
                    last_hit_idx: i,
                    direction: d,
                    gaps: 0,
                });
            }
            (Some(mut cs), Some(wd)) if cs.direction == wd => {
                let gap_windows = (i - cs.last_hit_idx).saturating_sub(1) as u32;
                cs.gaps += gap_windows;
                cs.last_hit_idx = i;
                current = Some(cs);
            }
            (Some(cs), Some(wd)) => {
                completed.push(cs);
                current = Some(RunState {
                    first_idx: i,
                    last_hit_idx: i,
                    direction: wd,
                    gaps: 0,
                });
            }
            (Some(cs), None) => {
                if i - cs.last_hit_idx > max_gap_windows {
                    completed.push(cs);
                    current = None;
                }
            }
        }
    }
    if let Some(cs) = current {
        completed.push(cs);
    }

    let mut runs: Vec<BoundaryRun> = completed
        .into_iter()
        .map(|rs| {
            let start = windows[rs.first_idx].start;
            let end = windows[rs.last_hit_idx].end;
            let length = end - start + 1;
            let slice = &windows[rs.first_idx..=rs.last_hit_idx];
            let mean_normalized_depth =
                slice.iter().map(|w| w.normalized_depth).sum::<f64>() / slice.len() as f64;
            BoundaryRun {
                start,
                end,
                length,
                mean_normalized_depth,
                direction: rs.direction,
                bridged_gaps: rs.gaps,
            }
        })
        .collect();

    runs.sort_by_key(|run| std::cmp::Reverse(run.length));
    runs
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mock::MockAlignmentSource;

    /// Test config with a small 300 kb boundary region for fast tests.
    /// The other regions from Default remain but are unused by boundary_scan.
    fn small_boundary_config() -> AnalysisConfig {
        AnalysisConfig {
            boundary_scan_region: Region::new("chr17", 1, 300_000),
            ..Default::default()
        }
    }

    #[test]
    fn boundary_scan_cn2_uniform_has_no_run() {
        let config = small_boundary_config();
        let mut mock = MockAlignmentSource::with_coverage(&config.boundary_scan_region, 30);
        let result = boundary_scan(&mut mock, 30.0, &config).unwrap();
        assert_eq!(result.windows.len(), 30);
        assert!(
            result.runs.is_empty(),
            "expected no run, got {:?}",
            result.longest_run()
        );
    }

    #[test]
    fn boundary_scan_cn3_cmt1a_detects_duplication() {
        let config = small_boundary_config();
        // 30× baseline across the whole region, +15× stacked on the middle
        // 100 kb = 45× effective depth there.
        let elevated = Region::new("chr17", 100_001, 200_000);
        let mut mock = MockAlignmentSource::with_coverage(&config.boundary_scan_region, 30)
            .add_coverage(&elevated, 15);
        let result = boundary_scan(&mut mock, 30.0, &config).unwrap();
        let run = result.longest_run().expect("expected a duplication run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        let tolerance: u32 = 30_000;
        assert!(
            run.length.abs_diff(100_000) <= tolerance,
            "run length {} not within 100_000 ± {}",
            run.length,
            tolerance
        );
        assert!(
            run.mean_normalized_depth > 1.3,
            "mean normalised depth {} should exceed 1.3",
            run.mean_normalized_depth
        );
    }

    #[test]
    fn boundary_scan_cn1_hnpp_detects_deletion() {
        let config = small_boundary_config();
        // Split the region into three flanks: left and right at 30×, middle
        // at 15× (half depth, simulating a deletion).
        let left = Region::new("chr17", 1, 100_000);
        let middle = Region::new("chr17", 100_001, 200_000);
        let right = Region::new("chr17", 200_001, 300_000);
        let mut mock = MockAlignmentSource::with_coverage(&left, 30)
            .add_coverage(&middle, 15)
            .add_coverage(&right, 30);
        let result = boundary_scan(&mut mock, 30.0, &config).unwrap();
        let run = result.longest_run().expect("expected a deletion run");
        assert_eq!(run.direction, BoundaryDirection::Deletion);
        let tolerance: u32 = 30_000;
        assert!(
            run.length.abs_diff(100_000) <= tolerance,
            "run length {} not within 100_000 ± {}",
            run.length,
            tolerance
        );
        assert!(
            run.mean_normalized_depth < 0.7,
            "mean normalised depth {} should be below 0.7",
            run.mean_normalized_depth
        );
    }

    #[test]
    fn boundary_scan_short_spike_is_filtered_by_min_length() {
        let config = small_boundary_config();
        // A single 10 kb window of elevated coverage. Below the default
        // min_reportable_run_length (100 kb) so it must be discarded.
        let spike = Region::new("chr17", 150_001, 160_000);
        let mut mock = MockAlignmentSource::with_coverage(&config.boundary_scan_region, 30)
            .add_coverage(&spike, 15);
        let result = boundary_scan(&mut mock, 30.0, &config).unwrap();
        assert!(
            result.runs.is_empty(),
            "10 kb spike should be below min_reportable_run_length"
        );
    }

    #[test]
    fn boundary_scan_makes_exactly_one_query_region_call() {
        let config = small_boundary_config();
        let mut mock = MockAlignmentSource::with_coverage(&config.boundary_scan_region, 30);
        boundary_scan(&mut mock, 30.0, &config).unwrap();
        assert_eq!(
            mock.query_count(),
            1,
            "boundary_scan must walk the region in a single query_region call"
        );
    }

    /// Helper to build a `WindowDepth` at a given index using 10 kb windows
    /// starting at position 1 with the specified normalised depth.
    fn window_at(index: u32, normalized: f64) -> WindowDepth {
        let start = 1 + index * WINDOW_SIZE;
        let end = start + WINDOW_SIZE - 1;
        WindowDepth {
            start,
            end,
            raw_depth: normalized * 30.0, // arbitrary; raw_depth isn't asserted on
            normalized_depth: normalized,
            unique_mapper_fraction: 1.0, // test windows default to fully unique
        }
    }

    #[test]
    fn find_longest_run_handles_direction_change() {
        // [Dup, Dup, Dup, Del, Del] - direction changes at index 3 without
        // passing through a Normal window first. The duplication run (length
        // 3) should win the longest-run comparison against the deletion run
        // (length 2). This exercises the (Some(d1), Some(d2)) state-machine
        // transition that the higher-level boundary_scan tests don't cover.
        let windows = vec![
            window_at(0, 1.5),
            window_at(1, 1.5),
            window_at(2, 1.5),
            window_at(3, 0.5),
            window_at(4, 0.5),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        assert_eq!(run.length, 30_000);
        assert_eq!(run.start, 1);
        assert_eq!(run.end, 30_000);
    }

    #[test]
    fn find_longest_run_picks_longer_direction_in_mixed_window_list() {
        // [Del, Del, Dup, Dup, Dup] - deletion comes first but is shorter.
        // Verifies the state machine doesn't depend on which direction the
        // window list happens to encounter first.
        let windows = vec![
            window_at(0, 0.5),
            window_at(1, 0.5),
            window_at(2, 1.5),
            window_at(3, 1.5),
            window_at(4, 1.5),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        assert_eq!(run.length, 30_000);
        assert_eq!(run.start, 20_001);
        assert_eq!(run.end, 50_000);
    }

    #[test]
    fn find_longest_run_earlier_wins_equal_length_tie() {
        // [Dup, Dup, Del, Del] - both runs have length 2. The earlier
        // (duplication) run should win per the documented tie-breaking rule.
        let windows = vec![
            window_at(0, 1.5),
            window_at(1, 1.5),
            window_at(2, 0.5),
            window_at(3, 0.5),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        assert_eq!(run.length, 20_000);
        assert_eq!(run.start, 1);
        assert_eq!(run.end, 20_000);
    }

    #[test]
    fn find_longest_run_bridges_small_gap() {
        // [Dup×3, Normal×3, Dup×3] - a 3-window gap is inside the bridge
        // budget (MAX_RUN_GAP_WINDOWS = 5), so the whole 90 kb span becomes
        // a single run. The mean normalised depth is diluted by the three
        // normal windows inside the span, which honestly reflects the
        // average depth across the CNV event.
        let windows = vec![
            window_at(0, 1.5),
            window_at(1, 1.5),
            window_at(2, 1.5),
            window_at(3, 1.0),
            window_at(4, 1.0),
            window_at(5, 1.0),
            window_at(6, 1.5),
            window_at(7, 1.5),
            window_at(8, 1.5),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a bridged run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        assert_eq!(run.start, 1);
        assert_eq!(run.end, 90_000);
        assert_eq!(run.length, 90_000);
        assert_eq!(run.bridged_gaps, 3, "should count 3 bridged gap windows");
        // Mean over nine windows: 6 × 1.5 + 3 × 1.0 = 12.0 / 9 = 1.333...
        assert!(
            (run.mean_normalized_depth - 1.333).abs() < 0.01,
            "mean {} should be ~1.333 (diluted by gap windows)",
            run.mean_normalized_depth
        );
    }

    #[test]
    fn find_longest_run_does_not_bridge_large_gap() {
        // [Dup×3, Normal×6, Dup×3] - a 6-window gap exceeds the budget of 5,
        // so the two segments are independent runs. Longer tie-break (both
        // length 30 kb) picks the earlier one.
        let windows = vec![
            window_at(0, 1.5),
            window_at(1, 1.5),
            window_at(2, 1.5),
            window_at(3, 1.0),
            window_at(4, 1.0),
            window_at(5, 1.0),
            window_at(6, 1.0),
            window_at(7, 1.0),
            window_at(8, 1.0),
            window_at(9, 1.5),
            window_at(10, 1.5),
            window_at(11, 1.5),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a run");
        assert_eq!(run.direction, BoundaryDirection::Duplication);
        assert_eq!(run.start, 1);
        assert_eq!(run.end, 30_000);
        assert_eq!(run.length, 30_000);
    }

    #[test]
    fn find_longest_run_trims_trailing_gap() {
        // [Dup×3, Normal×4] - the tail-end normal windows sit inside the
        // bridge budget but never reach another duplication. Those trailing
        // normals must NOT be counted as part of the run: length is 30 kb,
        // not 70 kb.
        let windows = vec![
            window_at(0, 1.5),
            window_at(1, 1.5),
            window_at(2, 1.5),
            window_at(3, 1.0),
            window_at(4, 1.0),
            window_at(5, 1.0),
            window_at(6, 1.0),
        ];
        let run = find_all_runs(&windows, 5, 1.3, 0.7, 0.0)
            .into_iter()
            .next()
            .expect("expected a run");
        assert_eq!(run.length, 30_000);
        assert_eq!(run.end, 30_000);
    }

    #[test]
    fn low_mappability_windows_excluded_from_run() {
        // Three "deletion" windows but all have low MAPQ - should NOT
        // produce a deletion call.
        let windows: Vec<WindowDepth> = (0..5)
            .map(|i| {
                let mut w = window_at(i, if (1..=3).contains(&i) { 0.5 } else { 1.0 });
                if (1..=3).contains(&i) {
                    w.unique_mapper_fraction = 0.1; // 10% unique = paralog zone
                }
                w
            })
            .collect();

        // With mappability filter at 0.5, all three deletion windows are
        // excluded, so no run should be found.
        let result = find_all_runs(&windows, 5, 1.3, 0.7, 0.5);
        assert!(
            result.is_empty(),
            "low-mappability deletion windows should not produce a run"
        );

        // Without filter (0.0), we do get a run (confirming the test setup works).
        let result_no_filter = find_all_runs(&windows, 5, 1.3, 0.7, 0.0);
        assert!(
            !result_no_filter.is_empty(),
            "without filter, run should exist"
        );
    }

    #[test]
    fn short_run_below_min_reportable_length_is_discarded() {
        // A 20 kb elevated region should NOT be reported as a duplication.
        // The 19.5 kb SMS-REP middle artifact must not produce a
        // classification.
        let mut config = small_boundary_config();
        config.min_reportable_run_length = 100_000; // 100 kb minimum
        let left = Region::new("chr17", 1, 100_000);
        let dup = Region::new("chr17", 100_001, 120_000); // only 20 kb elevated
        let right = Region::new("chr17", 120_001, 300_000);
        let mut mock = MockAlignmentSource::with_coverage(&left, 30)
            .add_coverage(&dup, 45) // 1.5x = duplication
            .add_coverage(&right, 30);
        let result = boundary_scan(&mut mock, 30.0, &config).unwrap();
        assert!(
            result.runs.is_empty(),
            "20 kb run should be below min_reportable_run_length (100 kb)"
        );
    }

    #[test]
    fn boundary_scan_rejects_non_positive_baseline() {
        let config = small_boundary_config();
        let mut mock = MockAlignmentSource::with_coverage(&config.boundary_scan_region, 30);
        let err = match boundary_scan(&mut mock, 0.0, &config) {
            Ok(r) => panic!("expected Internal error, got {r:?}"),
            Err(e) => e,
        };
        assert!(matches!(err, Error::Internal(_)));
        assert_eq!(err.exit_code(), 99);
    }
}
