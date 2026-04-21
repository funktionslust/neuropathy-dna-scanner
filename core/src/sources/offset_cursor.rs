//! `OffsetCursor` - a `Read + Seek` implementation over a sparse in-memory
//! buffer where real bytes live at known file offsets.
//!
//! Used by the WASM path: JS fetches specific byte ranges from a large
//! BAM/CRAM via `File.slice()`, then passes those ranges (with their
//! original file offsets) into WASM. noodles reads from the cursor as if
//! it were a contiguous file; seeks and reads within the prefetched ranges
//! work transparently, while reads into gaps return `UnexpectedEof`.

use std::io::{self, ErrorKind, Read, Seek, SeekFrom};

/// A sparse in-memory cursor that implements `Read + Seek` over
/// non-contiguous byte ranges placed at known file offsets.
///
/// Construct with `OffsetCursor::new(ranges)` where each range is
/// `(file_offset, bytes)`. Ranges are sorted on construction; overlapping
/// ranges cause a panic in debug builds.
pub struct OffsetCursor {
    /// Sorted by offset, non-overlapping.
    ranges: Vec<(u64, Vec<u8>)>,
    /// Optional: single backing buffer for the concatenated mode.
    /// When set, `ranges` entries hold empty vecs and data is read
    /// from `backing` using `range_offsets` to locate each chunk.
    backing: Option<Vec<u8>>,
    /// (start_in_backing, length) for each range when using backing buffer.
    range_offsets: Vec<(usize, usize)>,
    /// Current read/seek position in the virtual file.
    position: u64,
    /// Virtual file size: end of the last range.
    size: u64,
}

impl OffsetCursor {
    /// Construct a new cursor from a set of `(offset, bytes)` ranges.
    /// Ranges are sorted by offset. Overlapping ranges panic in debug.
    pub fn new(mut ranges: Vec<(u64, Vec<u8>)>) -> Self {
        ranges.sort_by_key(|(offset, _)| *offset);

        let range_offsets: Vec<(usize, usize)> = Vec::new();
        let size = ranges
            .last()
            .map(|(off, bytes)| off + bytes.len() as u64)
            .unwrap_or(0);
        Self {
            ranges,
            backing: None,
            range_offsets,
            position: 0,
            size,
        }
    }

    /// Construct from a single concatenated buffer + an offset table.
    /// Avoids copying: the `data` buffer is held as-is, and each range
    /// is a view into it. This is the zero-copy path for WASM where
    /// File.slice() results are concatenated in JS.
    pub fn from_concatenated(
        data: Vec<u8>,
        range_descs: &[(u64, u64)], // (file_offset, length)
    ) -> Self {
        let mut file_offsets: Vec<u64> = Vec::with_capacity(range_descs.len());
        let mut range_offsets: Vec<(usize, usize)> = Vec::with_capacity(range_descs.len());
        let mut buf_offset = 0usize;
        for &(file_offset, length) in range_descs {
            file_offsets.push(file_offset);
            range_offsets.push((buf_offset, length as usize));
            buf_offset += length as usize;
        }
        // Build ranges vec with file offsets + empty vecs. The actual data
        // is read from `backing` via `range_offsets`. The Vec<u8> in each
        // range is empty - only the file offset and length (via range_offsets)
        // are used for lookup.
        let ranges: Vec<(u64, Vec<u8>)> = file_offsets.iter().map(|&fo| (fo, Vec::new())).collect();
        let size = range_descs.last().map(|&(off, len)| off + len).unwrap_or(0);
        Self {
            ranges,
            backing: Some(data),
            range_offsets,
            position: 0,
            size,
        }
    }

    /// Find the range index covering `pos`, or `None` if `pos` is in a gap.
    fn find_range(&self, pos: u64) -> Option<usize> {
        let idx = self.ranges.partition_point(|(off, _)| *off <= pos);
        if idx == 0 {
            return None;
        }
        let i = idx - 1;
        let (off, _) = &self.ranges[i];
        let range_len = self.range_len(i);
        if pos < off + range_len as u64 {
            Some(i)
        } else {
            None
        }
    }

    /// Get the length of range `i`.
    fn range_len(&self, i: usize) -> usize {
        if self.backing.is_some() && i < self.range_offsets.len() {
            self.range_offsets[i].1
        } else {
            self.ranges[i].1.len()
        }
    }

    /// Get a byte slice for range `i`.
    fn range_bytes(&self, i: usize) -> &[u8] {
        if let Some(ref backing) = self.backing {
            if i < self.range_offsets.len() {
                let (start, len) = self.range_offsets[i];
                return &backing[start..start + len];
            }
        }
        &self.ranges[i].1
    }
}

impl Read for OffsetCursor {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let Some(i) = self.find_range(self.position) else {
            return Err(io::Error::new(
                ErrorKind::UnexpectedEof,
                format!(
                    "OffsetCursor: position {} is in a gap (no prefetched data)",
                    self.position
                ),
            ));
        };
        let (off, _) = &self.ranges[i];
        let bytes = self.range_bytes(i);
        let offset_in_range = (self.position - off) as usize;
        let available = bytes.len() - offset_in_range;
        let to_copy = buf.len().min(available);
        buf[..to_copy].copy_from_slice(&bytes[offset_in_range..offset_in_range + to_copy]);
        self.position += to_copy as u64;
        Ok(to_copy)
    }
}

impl Seek for OffsetCursor {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let new_pos = match pos {
            SeekFrom::Start(n) => n as i64,
            SeekFrom::Current(n) => self.position as i64 + n,
            SeekFrom::End(n) => self.size as i64 + n,
        };
        if new_pos < 0 {
            return Err(io::Error::new(
                ErrorKind::InvalidInput,
                "OffsetCursor: seek to negative position",
            ));
        }
        self.position = new_pos as u64;
        Ok(self.position)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_ranges() -> OffsetCursor {
        // Range A: offset 100, 10 bytes [0..9]
        // Range B: offset 200, 10 bytes [10..19]
        // Gap: 110..199
        OffsetCursor::new(vec![
            (100, (0u8..10).collect()),
            (200, (10u8..20).collect()),
        ])
    }

    #[test]
    fn read_within_single_range() {
        let mut c = two_ranges();
        c.seek(SeekFrom::Start(100)).unwrap();
        let mut buf = [0u8; 5];
        let n = c.read(&mut buf).unwrap();
        assert_eq!(n, 5);
        assert_eq!(buf, [0, 1, 2, 3, 4]);
    }

    #[test]
    fn read_across_range_boundary_gives_short_read() {
        let mut c = two_ranges();
        c.seek(SeekFrom::Start(105)).unwrap();
        let mut buf = [0u8; 20]; // asks for 20, only 5 available in range A
        let n = c.read(&mut buf).unwrap();
        assert_eq!(n, 5);
        assert_eq!(&buf[..5], &[5, 6, 7, 8, 9]);
    }

    #[test]
    fn read_in_gap_returns_unexpected_eof() {
        let mut c = two_ranges();
        c.seek(SeekFrom::Start(150)).unwrap();
        let mut buf = [0u8; 5];
        let err = c.read(&mut buf).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::UnexpectedEof);
    }

    #[test]
    fn seek_start_current_end() {
        let mut c = two_ranges();
        assert_eq!(c.seek(SeekFrom::Start(100)).unwrap(), 100);
        assert_eq!(c.seek(SeekFrom::Current(5)).unwrap(), 105);
        assert_eq!(c.seek(SeekFrom::Current(-3)).unwrap(), 102);
        // End is 200 + 10 = 210
        assert_eq!(c.seek(SeekFrom::End(0)).unwrap(), 210);
        assert_eq!(c.seek(SeekFrom::End(-10)).unwrap(), 200);
    }

    #[test]
    fn seek_into_gap_then_read_errors() {
        let mut c = two_ranges();
        c.seek(SeekFrom::Start(110)).unwrap(); // gap
        let mut buf = [0u8; 1];
        assert_eq!(
            c.read(&mut buf).unwrap_err().kind(),
            ErrorKind::UnexpectedEof
        );
    }

    #[test]
    fn read_second_range() {
        let mut c = two_ranges();
        c.seek(SeekFrom::Start(200)).unwrap();
        let mut buf = [0u8; 10];
        let n = c.read(&mut buf).unwrap();
        assert_eq!(n, 10);
        assert_eq!(buf, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]);
    }

    #[test]
    fn contiguous_single_range() {
        // A single range covering the entire "file" - no gaps.
        let data = vec![42u8; 100];
        let mut c = OffsetCursor::new(vec![(0, data)]);
        let mut buf = [0u8; 100];
        let n = c.read(&mut buf).unwrap();
        assert_eq!(n, 100);
        assert!(buf.iter().all(|&b| b == 42));
        assert_eq!(c.seek(SeekFrom::End(0)).unwrap(), 100);
    }

    #[test]
    fn empty_cursor() {
        let mut c = OffsetCursor::new(vec![]);
        let mut buf = [0u8; 1];
        let err = c.read(&mut buf).unwrap_err();
        assert_eq!(err.kind(), ErrorKind::UnexpectedEof);
        assert_eq!(c.seek(SeekFrom::End(0)).unwrap(), 0);
    }
}
