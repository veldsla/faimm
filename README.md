# faimm
Random access to indexed fasta using a memory mapped file.

# Usage
This crate provides indexed fasta access by using a memory mapped file to read
the sequence data. It is intended for accessing sequence data on genome sized
fasta files and provides random access based on base coordinates. Because an
indexed fasta file uses a limited number of bases per line separated by
(sometimes platform-specific) newlines you cannot directly use the bytes
available from the mmap.

Access is provided using a view of the mmap using zero-based base coordinates.
This view can then be used to iterate over bases (represented as `u8`) or
parsed into a string. Naive GC counting is also available.

Access to the sequence data doesn't require the `IndexedFasta` to be mutable.
This makes it easy to share.

## Example
```rust
use faimm::IndexedFasta;
let fa = IndexedFasta::from_file("test/genome.fa").expect("Error opening fa");
let chr_index = fa.fai().tid("ACGT-25").expect("Cannot find chr in index");
let v = fa.view(chr_index,0,50).expect("Cannot get .fa view");
//count the bases
let counts = v.count_bases();
//or print the sequence
println!("{}", v.to_string());
```

## Limitations
The parser uses a simple ASCII mask for allowable characters (64..128), does
not apply any IUPAC conversion or validation. Anything outside this range is
silently skipped. This means that also invalid `fasta` will be parsed. The mere
presence of an accompanying `.fai` provides the assumption of a valid fasta.
Requires Rust >=1.32

# Alternatives
[Rust-bio](https://crates.io/crates/bio) provides a competent indexed fasta
reader. The major difference is that it has an internal buffer an therefore
needs to be mutable when performing read operations. faimm is also faster. If
you want record based access (without an .fai index file)
[rust-bio](https://crates.io/crates/bio) or
[seq_io](https://crates.io/crates/seq_io) provide this.

# Performance
Calculating the GC content of target regions of an exome (231_410 regions) on
the Human reference (GRCh38) takes about 0.7 seconds (warm cache), slightly
faster than bedtools nuc (0.9s probably a more sound implementation) and
rust-bio (1.3s same implementation as example) Some tests show counting can
also be improved using SIMD, but nothing has been released.


