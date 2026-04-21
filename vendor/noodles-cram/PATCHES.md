# noodles-cram 0.92.0 vendor patches

This is a patched copy of [noodles-cram 0.92.0](https://crates.io/crates/noodles-cram/0.92.0).
Remove this patch once upstream publishes fixes.

## Changed files (7 files diffed against crates.io 0.92.0)

### 1. Huffman byte codec (`container/compression_header/encoding/codec/byte.rs`)
**Problem:** `decode_take` for Huffman-encoded bytes was `todo!()`, causing panics
on Sentieon-produced CRAMs that use this codec.

**Fix:** Implements single-symbol and multi-symbol Huffman byte decoding using
`CanonicalHuffmanDecoder`.

### 2. Huffman decoder (new file: `huffman.rs`)
New file implementing a canonical Huffman decoder used by the byte codec above.
Module declaration added to `lib.rs`.

### 3. Missing reference sequence handling (`io/reader/container/slice.rs`)
**Problem:** `get_record_reference_sequence()` panicked with
`expect("invalid reference sequence name")` when a chromosome's reference
wasn't loaded.

**Fix:** Returns `Ok(None)` instead of panicking. Records without reference
get `reference_sequence = None` -- positions are still available, only sequence
decoding is skipped.

### 4. Sequence iterator no-reference fix (`record/sequence/iter.rs`)
**Problem:** When iterating bases without a reference sequence, the code
returned `Some(b'N')` without advancing `last_read_position`, causing an
infinite loop and OOM.

**Fix:** Advances `last_read_position` by one per N emission, correctly
terminating after `read_length` bases.

### 5. Query skip-before optimization (`io/reader/query.rs`)
**Problem:** Query iterator read and decompressed every container from
chromosome start to the target region, even for queries at high positions.

**Fix:** Skips containers that end before the query range without seeking to
them. Adds early termination when containers start past the query end.

### 6. External data readers (`io/reader/container/slice/records/external_data_readers.rs`)
Added `alloc_owned_buf` method to support Huffman-decoded byte arrays that
need owned storage with a borrowed reference returned to the caller.

### 7. Module declaration (`lib.rs`)
Added `mod huffman;` to register the new Huffman decoder module.
