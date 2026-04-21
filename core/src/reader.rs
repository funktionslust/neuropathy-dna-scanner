//! Alignment I/O abstraction.
//!
//! [`AlignmentSource`] trait and supporting value types. Concrete
//! implementations:
//!
//! - `FileAlignmentSource` backed by `noodles::io::IndexedReader` over
//!   `std::fs::File`.
//! - `CursorAlignmentSource` backed by a sparse `Cursor<Vec<u8>>` populated
//!   from browser `File.slice()` range reads.
//! - [`crate::mock::MockAlignmentSource`] for unit tests.

use crate::Result;

/// Genomic region. Coordinates are **1-based inclusive** on both sides,
/// matching noodles' `Position` type.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
}

impl Region {
    /// Construct a region. Debug-asserts `start <= end` so inverted regions
    /// fail fast during development; release builds skip the check for zero
    /// overhead. 1-based inclusive means the smallest valid region has
    /// `start == end` and `length() == 1`.
    pub fn new(chrom: impl Into<String>, start: u32, end: u32) -> Self {
        debug_assert!(
            start <= end,
            "Region::new: start ({start}) must be <= end ({end}) - 1-based inclusive"
        );
        Self {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    /// Length of the region in bases. 1-based inclusive means
    /// `end - start + 1`.
    pub fn length(&self) -> u32 {
        self.end - self.start + 1
    }
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

/// Minimal alignment record. The trait layer abstracts away whether the
/// underlying format is CRAM or BAM; implementations must filter out
/// unmapped, secondary, and supplementary records before yielding them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentRecord {
    /// 1-based inclusive start position on the reference.
    pub start: u32,
    /// 1-based inclusive end position on the reference.
    pub end: u32,
    pub mapq: u8,
    /// `true` if this is a primary alignment. v0.1 implementations always set
    /// this to `true` because they filter non-primary records at source; the
    /// field exists so future implementations can optionally preserve
    /// secondary records for diagnostics without breaking the trait shape.
    pub is_primary: bool,
}

/// One `@SQ` line's worth of header metadata. `md5` comes from the `M5` tag
/// if the header provides it - it is used for chr17 MD5 verification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceInfo {
    pub name: String,
    pub length: u64,
    pub md5: Option<[u8; 16]>,
}

/// Parsed alignment header. v0.1 only carries the `@SQ` sequence list; other
/// header fields (program records, read groups, etc.) can be added later if
/// downstream stories need them.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AlignmentHeader {
    pub sequences: Vec<SequenceInfo>,
}

impl AlignmentHeader {
    /// Look up the chr17 MD5 from the `@SQ M5` tag. Accepts UCSC-style
    /// `"chr17"`, plain `"17"`, and any case variant (`Chr17`, `CHR17`, etc.)
    /// because BAM/CRAM headers are not case-normalised when written and
    /// real-world consumer WGS pipelines emit mixed case. Returns `None` if
    /// chr17 is absent from the header or if the present record has no `M5`
    /// tag.
    pub fn chr17_md5(&self) -> Option<[u8; 16]> {
        self.sequences
            .iter()
            .find(|s| s.name.eq_ignore_ascii_case("chr17") || s.name.eq_ignore_ascii_case("17"))
            .and_then(|s| s.md5)
    }
}

/// Source of alignment records for a given region. Implementations are
/// synchronous and I/O-capable: file-backed, cursor-backed for WASM, and
/// in-memory mock for unit tests
/// ([`crate::mock::MockAlignmentSource`]).
///
/// **Filter contract.** `query_region` must NOT yield unmapped, secondary, or
/// supplementary records. This filtering happens at the implementation layer
/// so the algorithm layer (depth_ratio, boundary_scan) can treat every record
/// it sees as a primary mapped alignment without additional checks.
///
/// **Overlap contract.** A record is "in" the region when its `[start, end]`
/// range intersects `region.[start, end]`, inclusive on both endpoints.
/// Concretely: `record.end >= region.start && record.start <= region.end`.
///
/// **Synchronous-only.** The WASM implementation prefetches bytes off-thread
/// in JS and then provides them synchronously to this trait via an
/// in-memory cursor.
pub trait AlignmentSource {
    /// Read and return the alignment file header. Safe to call multiple
    /// times; implementations may cache internally.
    fn read_header(&mut self) -> Result<AlignmentHeader>;

    /// Query records overlapping `region`. The returned iterator borrows from
    /// `self` so implementations can stream records without collecting into
    /// an intermediate `Vec`. Each yielded item is a `Result` because
    /// streaming implementations may fail part-way through a region (e.g.
    /// truncated CRAM container).
    #[allow(clippy::type_complexity)]
    fn query_region<'a>(
        &'a mut self,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<AlignmentRecord>> + 'a>>;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chr17_md5_is_case_insensitive() {
        let md5 = [0x12; 16];
        for name in ["chr17", "Chr17", "CHR17", "17"] {
            let header = AlignmentHeader {
                sequences: vec![SequenceInfo {
                    name: name.to_string(),
                    length: 83_257_441,
                    md5: Some(md5),
                }],
            };
            assert_eq!(
                header.chr17_md5(),
                Some(md5),
                "chr17_md5() failed for @SQ SN:{name}"
            );
        }
    }

    #[test]
    fn chr17_md5_returns_none_for_other_chromosomes() {
        let header = AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr1".to_string(),
                length: 248_956_422,
                md5: Some([0x34; 16]),
            }],
        };
        assert_eq!(header.chr17_md5(), None);
    }

    #[test]
    fn chr17_md5_returns_none_when_md5_tag_missing() {
        let header = AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: 83_257_441,
                md5: None,
            }],
        };
        assert_eq!(header.chr17_md5(), None);
    }

    #[test]
    fn region_length_is_inclusive() {
        assert_eq!(Region::new("chr17", 1000, 1000).length(), 1);
        assert_eq!(Region::new("chr17", 1000, 2000).length(), 1001);
    }

    #[test]
    fn region_display_is_ucsc_style() {
        let r = Region::new("chr17", 15_229_777, 15_265_079);
        assert_eq!(r.to_string(), "chr17:15229777-15265079");
        // Shorter form should still work.
        assert_eq!(Region::new("chr2", 1, 1).to_string(), "chr2:1-1");
    }
}
