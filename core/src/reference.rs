//! Reference fasta abstraction.
//!
//! [`ReferenceSource`] trait. Concrete implementations:
//!
//! - File-backed via `noodles::fasta::io::IndexedReader`.
//! - Cursor-backed over prefetched FASTA bytes for WASM.
//! - Real chr17 MD5 computation using the `md-5` crate.
//! - [`crate::mock::MockReferenceSource`] for unit tests.

use crate::reader::Region;
use crate::Result;

/// Source of reference sequence bytes plus the chr17 MD5 used for alignment
/// header verification. Both methods take `&mut self` because noodles'
/// indexed FASTA readers advance an internal cursor on every fetch.
pub trait ReferenceSource {
    /// Fetch the sequence bytes covering `region`. Output is uppercase ASCII
    /// IUPAC (`A`, `C`, `G`, `T`, `N`, plus ambiguity codes). The returned
    /// `Vec<u8>` has length `region.length()`.
    fn fetch_sequence(&mut self, region: &Region) -> Result<Vec<u8>>;

    /// Compute the MD5 of the full chr17 sequence for comparison against the
    /// `@SQ M5` tag in the alignment header.
    fn chr17_md5(&mut self) -> Result<[u8; 16]>;
}
