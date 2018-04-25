//! This crate provides an indexed fasta reader that uses a memory mapped file to read the sequence
//! data. It is intended for accessing sequence data on genome sized fasta files and provides
//! random access based on base coordinates. Because an indexed fasta file uses a limited number of
//! bases per line separated by (sometimes platform-specific) newlines you cannot directly use the
//! bytes available from the mmap.
//!
//! Access is provided using a view of the mmap using zero-based base coordinates. This view can
//! then be used to iterate over bases (represented as `u8`) or parsed into a string.
//!
//! Access to the sequence data doesn't require the `IndexedFastaReader` to be mutable. This makes
//! it easy to share.
//!
//! # Example
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
//! reference (GRCh38) takes about 0.7 seconds, slightly faster than bedtools (0.9s probably a more
//! sound implementation) and rust-bio (1.3s same implementation as example)
//!
extern crate memmap;

use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use memmap::{MmapOptions, Mmap};

/// The object that stores the parsed fasta index file. You can use it to map chromosome names to
/// indexes and lookup offsets for chr-start:end coordinates
#[derive(Debug, Clone)]
pub struct Fai {
    chromosomes: Vec<FaiRecord>,
    name_map: HashMap<String, usize>
}

impl Fai {
    /// Open a fasta index file from path `P`.
    pub fn from_file<P: AsRef<Path>>(path: P) ->io::Result<Self> {
        let f = File::open(path)?;
        let br = BufReader::new(f);

        let mut name_map = HashMap::new();
        let mut chromosomes = Vec::new();

        for l in br.lines() {
            let line = l?;
            let p: Vec<_> = line.split('\t').collect();
            //FIXME provide custom errors
            assert_eq!(p.len(), 5);

            name_map.insert(p[0].to_owned(), chromosomes.len());

            let ioerr = |e, msg| {
                io::Error::new(io::ErrorKind::InvalidData, format!("{}:{}", msg, e))
            };
            chromosomes.push( FaiRecord {
                len: p[1].parse().map_err(|e| ioerr(e, "Error parsing chr len in .fai"))?,
                offset: p[2].parse().map_err(|e| ioerr(e, "Error parsing chr offset in .fai"))?,
                line_bases: p[3].parse().map_err(|e| ioerr(e, "Error parsing chr line_bases in .fai"))?,
                line_width: p[4].parse().map_err(|e| ioerr(e, "Error parsing chr line_width in .fai"))?,
            });
        }

        Ok(Fai { chromosomes, name_map })
    }

    /// Calculate the slice coordinates (byte offsets).
    /// tid is the index of the chromosome (lookup with `Fai::tid` if necessary.
    /// start, end: zero based coordinates of the requested range.
    ///
    /// Returns an tuple (start, end) if successful. `io::Error` otherwise.
    #[inline]
    pub fn offset(&self, tid: usize, start: usize, stop: usize) -> io::Result<(usize, usize)> {
        let chr = &self.chromosomes.get(tid)
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds"))?;
        if stop > chr.len {
            return Err(io::Error::new(io::ErrorKind::Other, "FASTA read interval was out of bounds"));
        }

        let start_offset = chr.offset
           + (start / chr.line_bases) * chr.line_width
           + start % chr.line_bases;
        let stop_offset = chr.offset
           + (stop / chr.line_bases) * chr.line_width
           + stop % chr.line_bases;

        Ok((start_offset, stop_offset))
    }
    
    /// Calculate the slice coordinates (byte offsets).
    /// tid is the index of the chromosome (lookup with `Fai::tid` if necessary.
    ///
    /// Returns an tuple (start, end) if successful. `io::Error` otherwise.
    #[inline]
    pub fn offset_tid(&self, tid: usize) -> io::Result<(usize, usize)> {
        let chr = &self.chromosomes.get(tid)
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Chromomsome tid was out of bounds"))?;
        let start_offset = chr.offset;
        let stop_offset = chr.offset
           + (chr.len / chr.line_bases) * chr.line_width
           + chr.len % chr.line_bases;
        Ok((start_offset, stop_offset))
    }
    
    /// Return the index of the chromosome by name in the fasta index.
    ///
    /// Returns the position of chr `name` if succesful, None otherwise.
    #[inline]
    pub fn tid(&self, name: &str) -> Option<usize> {
        self.name_map.get(name).cloned()
    }

    /// Return the names chromosomes in the fasta index in no particular order use `Fai::tid` to
    /// map to index.
    ///
    /// Returns a `Vec<&str>` with the chromosome names.
    pub fn names(&self) -> Vec<&str> {
        self.name_map.keys().map(|s| s.as_ref()).collect()
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

/// The `IndexFastaReader` can be used to open a fasta file that has a valid .fai index file.
pub struct IndexedFastaReader {
    mmap: Mmap,
    fasta_index: Fai
}

impl IndexedFastaReader {
    /// Open a fasta file from path `P`. It is assumed that it has a valid .fai index file. The
    /// .fai file is created by appending .fai to the fasta file.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let mut fai_path = path.as_ref().as_os_str().to_owned();
        fai_path.push(".fai");
        let fasta_index = Fai::from_file(&fai_path)?;

        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(IndexedFastaReader { mmap, fasta_index })
    }

    /// Use tid, start and end to calculate a slice on the Fasta file. Use this view to iterate
    /// over the bases.
    ///
    /// Returns FastaView for the provided chromsome, start, end if successful, Error otherwise.
    pub fn view(&self, tid: usize, start: usize, stop: usize) -> io::Result<FastaView> {
        if start > stop {
            return Err(io::Error::new(io::ErrorKind::Other, "Invalid query interval"));
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
/// An iterator over the bases (as `u8`, so `b'C'`, `b'G'`, etc) in a `FastaView`.
pub struct Bases<'a>(std::slice::Iter<'a, u8>);

impl<'a> Iterator for Bases<'a> {
    type Item = u8;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(&b) = self.0.next() {
            if b == b'\r' || b == b'\n' {
                self.next()
            } else {
                Some(b)
            }
        } else {
            None
        }
    }
}

impl<'a> FastaView<'a> {
    /// Count the occurences of A, C, G, T, N, and other in the current view. This function does
    /// not differentiate between upper or lower case bases.
    ///
    /// Returns a `BasecCounts` object.
    pub fn count_bases(&self) -> BaseCounts {
        let mut bc: BaseCounts = Default::default();

        for b in self.0 {
            match *b {
                b'A' | b'a' => bc.a += 1,
                b'C' | b'c' => bc.c += 1,
                b'G' | b'g' => bc.g += 1,
                b'T' | b't' => bc.t += 1,
                b'N' | b'n' => bc.n += 1,
                b'\n' | b'\r' => {},
                _ => bc.other += 1,
            }
        }

        bc
    }

    /// Iterator over the bases in the current view. Bases are returned as `u8` representations of
    /// the `char`s in the fasta file.
    pub fn bases(&self) -> Bases {
        Bases(self.0.iter())
    }
}

/// Returns a newly allocated, utf-validated string with the sequence data in `Self`
impl<'a> ToString for FastaView<'a> {
    fn to_string(&self) -> String {
        String::from_utf8(self.bases().collect()).unwrap()
    }
}

/// Object that contains count occurences of the most common bases in DNA genome references: A, C, G,
/// T, N and other.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct BaseCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub other: usize
}

/// Initialize basecount with zeros
impl Default for BaseCounts {
    fn default() -> BaseCounts {
        BaseCounts {a: 0, c: 0, g: 0, t: 0, n: 0, other: 0}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open() {
        let ir = IndexedFastaReader::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.fai().names().len(), 3);
        assert_eq!(ir.fai().tid("ACGT-25"), Some(2));
        assert_eq!(ir.fai().tid("NotFound"), None);
    }

    #[test]
    fn view() {
        let ir = IndexedFastaReader::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.view(0, 0, 10).unwrap().to_string(), "AAAAAAAAAA");
        assert!(ir.view(0, 0, 11).is_err());

        assert_eq!(ir.view(2, 38, 62).unwrap().to_string(), "CCCCCCCCCCCCGGGGGGGGGGGG");
        assert!(ir.view(0, 120, 130).is_err());
    }

    #[test]
    fn view_tid() {
        let ir = IndexedFastaReader::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.view_tid(0).unwrap().to_string(),
            "AAAAAAAAAA");
        assert_eq!(ir.view_tid(1).unwrap().to_string(),
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        assert_eq!(ir.view_tid(2).unwrap().to_string(),
            "AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTT");
        assert!(ir.view_tid(3).is_err());
    }

    #[test]
    fn view_bases() {
        let ir = IndexedFastaReader::from_file("test/genome.fa").unwrap();
        let v = ir.view(2, 48, 52).unwrap();
        let mut b = v.bases();
        assert_eq!(b.next(), Some(b'C'));
        assert_eq!(b.next(), Some(b'C'));
        assert_eq!(b.next(), Some(b'G'));
        assert_eq!(b.next(), Some(b'G'));
        assert_eq!(b.next(), None);
    }

    #[test]
    fn view_counts() {
        let ir = IndexedFastaReader::from_file("test/genome.fa").unwrap();
        assert_eq!(ir.view(2, 48, 52).unwrap().count_bases(), BaseCounts {c:2,g:2, ..Default::default()});
    }
}
