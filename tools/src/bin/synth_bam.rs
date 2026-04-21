//! Synthetic BAM fixture generator.
//!
//! Produces deterministic BAM files at specified coverage profiles for
//! use in integration tests.
//!
//! Usage:
//!   synth_bam --output <path> --coverage chr17:15229777-15265079=30 --coverage chr2:50000000-70000000=30

use std::path::PathBuf;
use std::process;

use clap::Parser;
use nds_tools::synth::{self, CoverageSpec};

#[derive(Parser)]
#[command(
    name = "synth_bam",
    about = "Generate synthetic BAM fixtures at controlled coverage profiles"
)]
struct Cli {
    /// Output BAM file path. A .bai index is written alongside.
    #[arg(short, long)]
    output: PathBuf,

    /// Coverage specification in the format chr:start-end=depth.
    /// May be repeated for multiple regions.
    #[arg(short, long, required = true)]
    coverage: Vec<String>,

    /// PRNG seed for deterministic output. Same seed + same coverage
    /// arguments = byte-identical BAM.
    #[arg(short, long, default_value_t = 42)]
    seed: u64,
}

fn main() {
    let cli = Cli::parse();

    let specs: Vec<CoverageSpec> = cli
        .coverage
        .iter()
        .map(|s| {
            CoverageSpec::parse(s).unwrap_or_else(|e| {
                eprintln!("error: {e}");
                process::exit(1);
            })
        })
        .collect();

    eprintln!(
        "synth_bam: generating {} with {} coverage spec(s), seed={}",
        cli.output.display(),
        specs.len(),
        cli.seed
    );
    for spec in &specs {
        eprintln!(
            "  {}:{}-{} @ {}x",
            spec.chrom, spec.start, spec.end, spec.depth
        );
    }

    if let Err(e) = synth::generate_fixture(&cli.output, &specs, cli.seed) {
        eprintln!("error: {e}");
        process::exit(1);
    }

    eprintln!("synth_bam: done - {}", cli.output.display());
    let bai = cli.output.with_extension("bam.bai");
    if bai.exists() {
        eprintln!("synth_bam: index - {}", bai.display());
    }
}
