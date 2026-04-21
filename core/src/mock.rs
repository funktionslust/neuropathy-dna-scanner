//! Test-only mock implementations of [`AlignmentSource`] and
//! [`ReferenceSource`]. Gated by `#[cfg(test)]` at the `lib.rs` declaration
//! so this module exists only when building tests; `pub(crate)` visibility
//! makes the types reachable from every unit test in the crate via
//! `crate::mock::*`.
//!
//! The core algorithm unit tests build against these mocks.

use crate::reader::{AlignmentHeader, AlignmentRecord, AlignmentSource, Region};
use crate::reference::ReferenceSource;
use crate::Result;

/// In-memory alignment source. Filters by position only; the chromosome on
/// the queried region is ignored. The two regions used by `depth_ratio`
/// (chr17 PMP22 at ~15 Mb and chr2 control at 50–70 Mb) don't overlap in
/// position space, so records inserted for one are naturally excluded from
/// queries on the other.
pub struct MockAlignmentSource {
    header: AlignmentHeader,
    records: Vec<AlignmentRecord>,
    query_count: usize,
}

impl MockAlignmentSource {
    /// Construct from a canned header and a flat list of records.
    pub fn new(header: AlignmentHeader, records: Vec<AlignmentRecord>) -> Self {
        Self {
            header,
            records,
            query_count: 0,
        }
    }

    /// Attach or replace the header. Builder-style for test ergonomics.
    pub fn with_header(mut self, header: AlignmentHeader) -> Self {
        self.header = header;
        self
    }

    /// Number of times `query_region` has been called on this mock since
    /// construction. Used to assert single-pass
    /// contracts (e.g. `boundary_scan` must call `query_region` exactly once).
    pub fn query_count(&self) -> usize {
        self.query_count
    }

    /// Synthesise uniform coverage at `depth` × across `region`. Record
    /// count is `region.length() * depth / READ_LENGTH`; reads step by
    /// `READ_LENGTH / depth` so the mean depth across the region is roughly
    /// `depth`. The `.max(1)` guard covers the zero-step corner case for
    /// very short regions or very high depths.
    pub fn with_coverage(region: &Region, depth: u32) -> Self {
        const READ_LENGTH: u32 = 150;
        assert!(depth > 0, "depth must be positive");
        let region_length = region.length();
        let record_count = region_length.saturating_mul(depth) / READ_LENGTH;
        let wrap = region_length.saturating_sub(READ_LENGTH - 1).max(1);
        let step = (READ_LENGTH / depth).max(1);
        let mut records = Vec::with_capacity(record_count as usize);
        for i in 0..record_count {
            let offset = (i * step) % wrap;
            let start = region.start + offset;
            let end = start + READ_LENGTH - 1;
            records.push(AlignmentRecord {
                start,
                end,
                mapq: 60,
                is_primary: true,
            });
        }
        Self {
            header: AlignmentHeader::default(),
            records,
            query_count: 0,
        }
    }

    /// Append synthetic coverage for another region. Useful for multi-region
    /// test setups - e.g. `with_coverage(pmp22, 45).add_coverage(chr2_ctl, 30)`
    /// for a CN=3 sample in the `depth_ratio` tests.
    pub fn add_coverage(mut self, region: &Region, depth: u32) -> Self {
        let more = Self::with_coverage(region, depth);
        self.records.extend(more.records);
        self
    }
}

impl AlignmentSource for MockAlignmentSource {
    fn read_header(&mut self) -> Result<AlignmentHeader> {
        Ok(self.header.clone())
    }

    fn query_region<'a>(
        &'a mut self,
        region: &Region,
    ) -> Result<Box<dyn Iterator<Item = Result<AlignmentRecord>> + 'a>> {
        self.query_count += 1;
        let rstart = region.start;
        let rend = region.end;
        let iter = self
            .records
            .iter()
            .filter(|r| r.is_primary)
            .filter(move |r| r.end >= rstart && r.start <= rend)
            .copied()
            .map(Ok);
        Ok(Box::new(iter))
    }
}

/// In-memory reference source. Returns empty sequence bytes (unit tests
/// don't exercise sequence fetches) and a canned chr17 MD5.
pub struct MockReferenceSource {
    chr17_md5: [u8; 16],
}

impl MockReferenceSource {
    pub fn new(chr17_md5: [u8; 16]) -> Self {
        Self { chr17_md5 }
    }
}

impl ReferenceSource for MockReferenceSource {
    fn fetch_sequence(&mut self, _region: &Region) -> Result<Vec<u8>> {
        Ok(Vec::new())
    }

    fn chr17_md5(&mut self) -> Result<[u8; 16]> {
        Ok(self.chr17_md5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::SequenceInfo;

    #[test]
    fn mock_alignment_source_returns_header() {
        let header = AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: 83_257_441,
                md5: Some([0xAB; 16]),
            }],
        };
        let mut mock = MockAlignmentSource::new(header.clone(), Vec::new());
        assert_eq!(mock.read_header().unwrap(), header);
    }

    #[test]
    fn mock_alignment_source_filters_by_region() {
        let records = vec![
            AlignmentRecord {
                start: 1000,
                end: 1149,
                mapq: 60,
                is_primary: true,
            },
            AlignmentRecord {
                start: 2000,
                end: 2149,
                mapq: 60,
                is_primary: true,
            },
            AlignmentRecord {
                start: 3000,
                end: 3149,
                mapq: 60,
                is_primary: true,
            },
        ];
        let mut mock = MockAlignmentSource::new(AlignmentHeader::default(), records);
        let region = Region::new("chr17", 1500, 2500);
        let returned: Vec<_> = mock
            .query_region(&region)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            returned.len(),
            1,
            "expected exactly one record in 1500..=2500"
        );
        assert_eq!(returned[0].start, 2000);
    }

    #[test]
    fn mock_alignment_source_with_coverage_produces_expected_read_count() {
        let region = Region::new("chr17", 1, 10_000);
        let mut mock = MockAlignmentSource::with_coverage(&region, 30);
        let returned: Vec<_> = mock
            .query_region(&region)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        // 10_000 bases × 30× coverage / 150 bp reads = 2_000 records
        assert_eq!(returned.len(), 2_000);
    }

    #[test]
    fn mock_reference_source_returns_chr17_md5() {
        let expected = [0xCD; 16];
        let mut mock = MockReferenceSource::new(expected);
        assert_eq!(mock.chr17_md5().unwrap(), expected);
    }

    #[test]
    fn with_header_replaces_default_header() {
        let replacement = AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr1".to_string(),
                length: 248_956_422,
                md5: Some([0xEF; 16]),
            }],
        };
        let mut mock = MockAlignmentSource::with_coverage(&Region::new("chr17", 1, 1_000), 10)
            .with_header(replacement.clone());
        assert_eq!(mock.read_header().unwrap(), replacement);
    }

    #[test]
    fn mock_filters_out_non_primary_records() {
        // Enforces the AlignmentSource filter contract: non-primary records
        // must not be yielded.
        let records = vec![
            AlignmentRecord {
                start: 1000,
                end: 1149,
                mapq: 60,
                is_primary: true,
            },
            AlignmentRecord {
                start: 1000,
                end: 1149,
                mapq: 60,
                is_primary: false,
            },
        ];
        let mut mock = MockAlignmentSource::new(AlignmentHeader::default(), records);
        let returned: Vec<_> = mock
            .query_region(&Region::new("chr17", 500, 2000))
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(returned.len(), 1);
        assert!(returned[0].is_primary);
    }

    #[test]
    fn read_header_is_safe_to_call_multiple_times() {
        let header = AlignmentHeader {
            sequences: vec![SequenceInfo {
                name: "chr17".to_string(),
                length: 83_257_441,
                md5: Some([0xAB; 16]),
            }],
        };
        let mut mock = MockAlignmentSource::new(header.clone(), Vec::new());
        let first = mock.read_header().unwrap();
        let second = mock.read_header().unwrap();
        assert_eq!(first, second);
        assert_eq!(first, header);
    }

    #[test]
    fn query_count_increments_on_each_call() {
        let mut mock = MockAlignmentSource::new(AlignmentHeader::default(), Vec::new());
        assert_eq!(mock.query_count(), 0);
        let _ = mock.query_region(&Region::new("chr17", 1, 100)).unwrap();
        assert_eq!(mock.query_count(), 1);
        let _ = mock.query_region(&Region::new("chr2", 1, 100)).unwrap();
        assert_eq!(mock.query_count(), 2);
        let _ = mock.query_region(&Region::new("chr17", 1, 100)).unwrap();
        assert_eq!(mock.query_count(), 3);
    }

    #[test]
    fn add_coverage_combines_records_from_multiple_regions() {
        let pmp22 = Region::new("chr17", 15_229_777, 15_229_777 + 999);
        let chr2_ctl = Region::new("chr2", 50_000_000, 50_000_000 + 999);
        let mut mock = MockAlignmentSource::with_coverage(&pmp22, 30).add_coverage(&chr2_ctl, 30);

        let pmp22_records: Vec<_> = mock
            .query_region(&pmp22)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();
        let chr2_records: Vec<_> = mock
            .query_region(&chr2_ctl)
            .unwrap()
            .collect::<Result<Vec<_>>>()
            .unwrap();

        // Each region: 1000 bases × 30× / 150 bp = 200 records.
        assert_eq!(pmp22_records.len(), 200);
        assert_eq!(chr2_records.len(), 200);
    }
}
