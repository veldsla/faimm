//! This crate provides indexed fasta access by using a memory mapped file to read the sequence
//! data. It is intended for accessing sequence data on genome sized fasta files and provides
//! random access based on base coordinates. Because an indexed fasta file uses a limited number of
//! bases per line separated by (sometimes platform-specific) newlines you cannot directly use the
//! bytes available from the mmap.
//!
//! Access is provided using a view of the mmap using zero-based base coordinates. This view can
//! then be used to iterate over bases (represented as `u8`) or parsed into a string. Naive gc
//! counting is also available.
//!
//! Access to the sequence data doesn't require the `IndexedFasta` to be mutable. This makes
//! it easy to share.
//!
//! # Example
//! ```
//! use faimm::IndexedFasta;
//! let fa = IndexedFasta::from_file("test/genome.fa").expect("Error opening fa");
//! let chr_index = fa.fai().tid("ACGT-25").expect("Cannot find chr in index");
//! let v = fa.view(chr_index,0,50).expect("Cannot get .fa view");
//! //count the bases
//! let counts = v.count_bases();
//! //or print the sequence
//! println!("{}", v.to_string());
//! ```
//! # Limitations
//! The parser uses a simple ascii mask for allowable characters (64..128), does not apply any
//! IUPAC converson or validation. Anything outside this range is silently skipped. This means that
//! also invalid `fasta` will be parsed. The mere presence of an accompanying `.fai` provides the
//! assumption of a valid fasta.
//! Requires Rust >=1.64
//!
//! # Alternatives
//! [Rust-bio](https://crates.io/crates/bio) provides a competent indexed fasta reader. The major
//! difference is that it has an internal buffer an therefore needs to be mutable when performing
//! read operations. faimm is also faster. If you want record based access (without an .fai index
//! file) [rust-bio](https://crates.io/crates/bio) or [seq_io](https://crates.io/crates/seq_io)
//! provide this.
//!
//! # Performance
//! Calculating the gc content of target regions of an exome (231_410 regions) on the Human
//! reference (GRCh38) takes about 0.7 seconds (warm cache), slightly faster than bedtools nuc (0.9s probably a more
//! sound implementation) and rust-bio (1.3s same implementation as example)
//! Some tests show counting can also be improved using simd, but nothing has been released.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use indexmap::IndexSet;
use memmap2::{Mmap, MmapOptions};

/// The object that stores the parsed fasta index file. You can use it to map chromosome names to
/// indexes and lookup offsets for chr-start:end coordinates
#[derive(Debug, Clone)]
pub struct Fai {
    chromosomes: Vec<FaiRecord>,
    name_map: IndexSet<String>,
}

impl Fai {
    /// Open a fasta index file from path `P`.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let f = File::open(path)?;
        let br = BufReader::new(f);

        let mut name_map = IndexSet::new();
        let mut chromosomes = Vec::new();

        for l in br.lines() {
            let line = l?;
            let p: Vec<_> = line.split('\t').collect();

            if p.len() != 5 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Expected 5 columns in .fai file.",
                ));
            }

            name_map.insert(p[0].to_owned());

            let ioerr =
                |e, msg| io::Error::new(io::ErrorKind::InvalidData, format!("{}:{}", msg, e));
            chromosomes.push(FaiRecord {
                len: p[1]
                    .parse()
                    .map_err(|e| ioerr(e, "Error parsing chr len in .fai"))?,
                offset: p[2]
                    .parse()
                    .map_err(|e| ioerr(e, "Error parsing chr offset in .fai"))?,
                line_bases: p[3]
                    .parse()
                    .map_err(|e| ioerr(e, "Error parsing chr line_bases in .fai"))?,
                line_width: p[4]
                    .parse()
                    .map_err(|e| ioerr(e, "Error parsing chr line_width in .fai"))?,
            });
        }

        Ok(Fai {
            chromosomes,
            name_map,
        })
    }

    /// Calculate the slice coordinates (byte offsets).
    /// tid is the index of the chromosome (lookup with `Fai::tid` if necessary.
    /// start, end: zero based coordinates of the requested range.
    ///
    /// Returns an tuple (start, end) if successful. `io::Error` otherwise.
    #[inline]
    pub fn offset(&self, tid: usize, start: usize, stop: usize) -> io::Result<(usize, usize)> {
        let chr = &self.chromosomes.get(tid).ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds")
        })?;
        if stop > chr.len {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "FASTA read interval was out of bounds",
            ));
        }

        let start_offset =
            chr.offset + (start / chr.line_bases) * chr.line_width + start % chr.line_bases;
        let stop_offset =
            chr.offset + (stop / chr.line_bases) * chr.line_width + stop % chr.line_bases;

        Ok((start_offset, stop_offset))
    }

    /// Calculate the slice coordinates (byte offsets).
    /// tid is the index of the chromosome (lookup with `Fai::tid` if necessary.
    ///
    /// Returns an tuple (start, end) if successful. `io::Error` otherwise.
    #[inline]
    pub fn offset_tid(&self, tid: usize) -> io::Result<(usize, usize)> {
        let chr = &self.chromosomes.get(tid).ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds")
        })?;
        let start_offset = chr.offset;
        let stop_offset =
            chr.offset + (chr.len / chr.line_bases) * chr.line_width + chr.len % chr.line_bases;
        Ok((start_offset, stop_offset))
    }

    /// Return the index of the chromosome by name in the fasta index.
    ///
    /// Returns the position of chr `name` if succesful, None otherwise.
    #[inline]
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.name_map.get_index_of(name)
    }

    /// Return the index of a chromosome in the fasta index.
    ///
    /// Returns the size in bases as usize.
    pub fn size(&self, tid: usize) -> io::Result<usize> {
        let chr = &self.chromosomes.get(tid).ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds")
        })?;
        Ok(chr.len)
    }

    /// Return the name of the chromomsome at index tid
    pub fn name(&self, tid: usize) -> io::Result<&String> {
        self.name_map.get_index(tid).ok_or_else(|| {
            io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds")
        })
    }

    /// Return the names of the chromosomes from the fasta index in the same order as in the
    /// `.fai`. You can use `Fai::tid` to map it back to an index.
    ///
    /// Returns a `Vec<&str>` with the chromosome names.
    pub fn names(&self) -> Vec<&str> {
        self.name_map.iter().map(|s| s.as_str()).collect()
    }
}

/// FaiRecord stores the length, offset, and fasta file characterics of a single chromosome
#[derive(Debug, Clone)]
pub struct FaiRecord {
    len: usize,
    offset: usize,
    line_bases: usize,
    line_width: usize,
}

/// The `IndexFasta` can be used to open a fasta file that has a valid .fai index file.
pub struct IndexedFasta {
    mmap: Mmap,
    fasta_index: Fai,
}

impl IndexedFasta {
    /// Open a fasta file from path `P`. It is assumed that it has a valid .fai index file. The
    /// .fai file is created by appending .fai to the fasta file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut fai_path = path.as_ref().as_os_str().to_owned();
        fai_path.push(".fai");
        let fasta_index = Fai::from_file(&fai_path)?;

        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(IndexedFasta { mmap, fasta_index })
    }

    /// Use tid, start and end to calculate a slice on the Fasta file. Use this view to iterate
    /// over the bases.
    ///
    /// Returns FastaView for the provided chromsome, start, end if successful, Error otherwise.
    pub fn view(&self, tid: usize, start: usize, stop: usize) -> io::Result<FastaView> {
        if start > stop {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid query interval",
            ));
        }

        let (start_byte, stop_byte) = self.fasta_index.offset(tid, start, stop)?;
        //println!("offset for chr {}:{}-{} is {}-{}", tid, start, stop, start_byte, stop_byte);
        Ok(FastaView(&self.mmap[start_byte..stop_byte]))
    }

    /// Use tid to return a view of an entire chromosome.
    ///
    /// Returns FastaView for the provided chromsome indicated by tid if successful, Error otherwise.
    pub fn view_tid(&self, tid: usize) -> io::Result<FastaView> {
        let (start_byte, stop_byte) = self.fasta_index.offset_tid(tid)?;
        //println!("offset for chr {}:{}-{} is {}-{}", tid, start, stop, start_byte, stop_byte);
        Ok(FastaView(&self.mmap[start_byte..stop_byte]))
    }

    /// Return a reference to the `Fai` that contains information from the fasta index.
    ///
    /// Returns a reference to `Fai`.
    pub fn fai(&self) -> &Fai {
        &self.fasta_index
    }
}

/// A view of a slice of the fasta file bounded by provided coordinates
pub struct FastaView<'a>(&'a [u8]);

impl<'a> FastaView<'a> {
    /// Count the occurences of A, C, G, T, N, and other in the current view. This function does
    /// not differentiate between upper or lower case bases.
    ///
    /// Returns a `BasecCounts` object.
    pub fn count_bases(&self) -> BaseCounts {
        let mut bc: BaseCounts = Default::default();

        for b in self.bases() {
            let v: u8 = b << 3;
            if v ^ 8 == 0 {
                bc.a += 1;
            } else if v ^ 24 == 0 {
                bc.c += 1;
            } else if v ^ 56 == 0 {
                bc.g += 1;
            } else if v ^ 112 == 0 {
                bc.n += 1;
            } else if v ^ 160 == 0 {
                bc.t += 1;
            } else {
                bc.other += 1;
            }
        }

        bc
    }

    /// Iterator over the bases in the current view. Bases are returned as `u8` representations of
    /// the `char`s in the fasta file. Keep only that chars between 164 and 128 (effectively
    /// skipping newlines)
    pub fn bases(&self) -> impl Iterator<Item = &'a u8> {
        self.0.iter().filter(|&&b| b & 192 == 64)
    }
}

/// Returns a newly allocated, utf8-validated string with the sequence data in `Self`
impl<'a> ToString for FastaView<'a> {
    fn to_string(&self) -> String {
        String::from_utf8(self.bases().cloned().collect()).unwrap()
    }
}

impl<'a> Read for FastaView<'a> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut read = 0;
        let mut skipped = 0;
        for (t, s) in buf.iter_mut().zip(self.0.iter().filter(|&&c| {
            let base = c & 192 == 64;
            if !base {
                skipped += 1;
            }
            base
        })) {
            *t = *s;
            read += 1;
        }
        self.0 = &self.0[(skipped + read)..];
        Ok(read)
    }
}

/// Object that contains count occurrences of the most common bases in DNA genome references: A, C, G,
/// T, N and other.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub other: usize,
}

/// Initialize basecount with zeros
impl Default for BaseCounts {
    fn default() -> BaseCounts {
        BaseCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            other: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fai() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.fai().names().len(), 3);
        assert_eq!(ir.fai().tid("ACGT-25"), Some(2));
        assert_eq!(ir.fai().tid("NotFound"), None);

        assert_eq!(ir.fai().size(2).unwrap(), 100);
        assert_eq!(ir.fai().name(2).unwrap(), "ACGT-25");
        assert!(ir.fai().name(3).is_err());
    }

    #[test]
    fn view() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.view(0, 0, 10).unwrap().to_string(), "AAAAAAAAAA");
        assert!(ir.view(0, 0, 11).is_err());

        assert_eq!(
            ir.view(2, 38, 62).unwrap().to_string(),
            "CCCCCCCCCCCCGGGGGGGGGGGG"
        );
        assert_eq!(
            ir.view(2, 74, 100).unwrap().to_string(),
            "GTTTTTTTTTTTTTTTTTTTTTTTTT"
        );
        assert!(ir.view(0, 120, 130).is_err());
    }

    #[test]
    fn view_tid() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.view_tid(0).unwrap().to_string(), "AAAAAAAAAA");
        assert_eq!(ir.view_tid(1).unwrap().to_string(),
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(ir.view_tid(2).unwrap().to_string(),
            "AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTT");
        assert!(ir.view_tid(3).is_err());
    }

    #[test]
    fn view_bases() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        let v = ir.view(2, 48, 52).unwrap();
        let mut b = v.bases();
        assert_eq!(b.next(), Some(&b'C'));
        assert_eq!(b.next(), Some(&b'C'));
        assert_eq!(b.next(), Some(&b'G'));
        assert_eq!(b.next(), Some(&b'G'));
        assert_eq!(b.next(), None);
    }

    #[test]
    fn view_counts() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        assert_eq!(
            ir.view(2, 48, 52).unwrap().count_bases(),
            BaseCounts {
                c: 2,
                g: 2,
                ..Default::default()
            }
        );
    }

    #[test]
    fn read_view() {
        let ir = IndexedFasta::from_file("test/genome.fa").unwrap();
        let mut buf = vec![0; 25];
        let mut v = ir.view_tid(2).unwrap();
        println!("{}", v.to_string());
        assert_eq!(v.read(&mut buf).unwrap(), 25);
        assert_eq!(buf, vec![b'A'; 25]);
        assert_eq!(v.read(&mut buf).unwrap(), 25);
        assert_eq!(buf, vec![b'C'; 25]);
        assert_eq!(v.read(&mut buf).unwrap(), 25);
        assert_eq!(buf, vec![b'G'; 25]);

        let mut buf2 = vec![0; 10];
        assert_eq!(v.read(&mut buf2).unwrap(), 10);
        assert_eq!(buf2, vec![b'T'; 10]);
        assert_eq!(v.read(&mut buf2).unwrap(), 10);
        assert_eq!(buf2, vec![b'T'; 10]);
        assert_eq!(v.read(&mut buf2).unwrap(), 5);
        assert_eq!(&buf2[0..5], vec![b'T'; 5].as_slice());
    }
}
