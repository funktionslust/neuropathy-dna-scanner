//! Convenience task runner for neuropathy-dna-scanner development.
//!
//! Usage:
//!   cargo xtask synth-fixtures    Generate all test fixtures

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("synth-fixtures") => synth_fixtures(),
        Some(other) => {
            eprintln!("unknown subcommand: {other}");
            eprintln!("usage: cargo xtask synth-fixtures");
            std::process::exit(1);
        }
        None => {
            eprintln!("usage: cargo xtask synth-fixtures");
            std::process::exit(1);
        }
    }
}

fn synth_fixtures() {
    // Resolve workspace root from CARGO_MANIFEST_DIR (tools/) -> parent.
    let workspace_root = std::env::var("CARGO_MANIFEST_DIR")
        .map(std::path::PathBuf::from)
        .ok()
        .and_then(|p| p.parent().map(|p| p.to_path_buf()))
        .unwrap_or_else(|| std::path::PathBuf::from("."));
    let out_dir = workspace_root.join("core/tests/data");
    std::fs::create_dir_all(&out_dir).expect("create output dir");

    for (name, specs) in nds_tools::synth::fixture_specs() {
        let path = out_dir.join(format!("{name}.bam"));
        eprintln!("xtask: generating {}", path.display());

        if let Err(e) = nds_tools::synth::generate_fixture(&path, &specs, 42) {
            eprintln!("xtask: failed for {name}: {e}");
            std::process::exit(1);
        }
    }

    // Error-path fixtures.
    eprintln!("xtask: generating error-path fixtures...");
    if let Err(e) = nds_tools::synth::generate_error_fixtures(&out_dir) {
        eprintln!("xtask: error fixtures failed: {e}");
        std::process::exit(1);
    }

    eprintln!("xtask: all fixtures generated in {}", out_dir.display());
}
