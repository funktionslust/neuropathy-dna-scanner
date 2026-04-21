use std::collections::HashMap;

use crate::container::block;

pub struct ExternalDataReaders<'c> {
    low_readers: [Option<&'c [u8]>; 64],
    high_readers: HashMap<block::ContentId, &'c [u8]>,
    /// Owned buffers for Huffman-decoded byte data. When a Huffman byte
    /// codec is used in `decode_take`, the decoded bytes are allocated
    /// here. Returned slices point into each Vec's heap allocation, which
    /// remains stable even if the outer Vec reallocates.
    owned_bufs: Vec<Vec<u8>>,
}

impl<'c> ExternalDataReaders<'c> {
    pub fn new() -> Self {
        Self {
            low_readers: init_low_readers(),
            high_readers: HashMap::new(),
            owned_bufs: Vec::new(),
        }
    }

    /// Allocate a buffer for Huffman-decoded data and return a borrowed
    /// slice into it. The buffer lives as long as this `ExternalDataReaders`.
    /// The returned pointer is stable: Vec's heap data does not move when
    /// the outer `owned_bufs` Vec reallocates (only the stack-side Vec
    /// metadata moves; the heap allocation stays put).
    pub fn alloc_owned_buf(&mut self, data: Vec<u8>) -> &'c [u8] {
        self.owned_bufs.push(data);
        let buf = self.owned_bufs.last().unwrap();
        // SAFETY: buf.as_ptr() points to a heap allocation that remains
        // valid for the lifetime of self (and thus 'c). The inner Vec's
        // heap data is not moved by outer-Vec reallocation.
        unsafe { std::slice::from_raw_parts(buf.as_ptr(), buf.len()) }
    }

    pub fn insert(&mut self, id: block::ContentId, reader: &'c [u8]) {
        match id {
            i @ 0..=63 => {
                self.low_readers[i as usize] = Some(reader);
            }
            _ => {
                self.high_readers.insert(id, reader);
            }
        }
    }

    pub fn get_mut(&mut self, id: &block::ContentId) -> Option<&mut &'c [u8]> {
        match *id {
            i @ 0..=63 => self.low_readers[i as usize].as_mut(),
            _ => self.high_readers.get_mut(id),
        }
    }
}

fn init_low_readers<'c>() -> [Option<&'c [u8]>; 64] {
    [
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None, None, None, None, None, None, None, None, None, None, None, None,
        None, None, None, None,
    ]
}
