//! Decompress the committed `core/data/catalog.json.gz` into `OUT_DIR` so
//! `lib.rs` can embed it via `include_bytes!`. Keeps the git-tracked blob
//! small while letting the build be self-contained on a fresh clone.

use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;

use flate2::read::GzDecoder;

fn main() {
    let manifest_dir = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let source = manifest_dir.join("../core/data/catalog.json.gz");
    let out_dir = PathBuf::from(std::env::var_os("OUT_DIR").unwrap());
    let target = out_dir.join("catalog.json");

    println!("cargo:rerun-if-changed={}", source.display());

    let input = File::open(&source).unwrap_or_else(|e| panic!("open {}: {e}", source.display()));
    let mut decoder = GzDecoder::new(input);
    let mut buf = Vec::new();
    decoder
        .read_to_end(&mut buf)
        .unwrap_or_else(|e| panic!("decompress {}: {e}", source.display()));

    let mut out =
        File::create(&target).unwrap_or_else(|e| panic!("create {}: {e}", target.display()));
    out.write_all(&buf)
        .unwrap_or_else(|e| panic!("write {}: {e}", target.display()));
}
