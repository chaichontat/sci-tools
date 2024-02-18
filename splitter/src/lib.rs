#![feature(portable_simd)]
#![feature(ascii_char)]

use std::{collections::HashMap, hash::BuildHasherDefault, io, path::Path};
use twox_hash::XxHash64;
use needletail::FastxReader;
pub mod simdstuffs;

pub struct FastqSequence {
    pub reader: Box<dyn FastxReader>,
}

impl Iterator for FastqSequence {
    type Item = (Vec<u8>, Vec<u8>);

    fn next(&mut self) -> Option<Self::Item> {
        let record = self.reader.next()?.expect("invalid record");
        Some((record.raw_seq().to_vec(), record.qual().unwrap().to_vec()))
    }
}

pub fn extract_idx(file: &Path, key_side: usize) -> io::Result<HashMap<String, String, BuildHasherDefault<XxHash64>>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(file)?;
    let mut map: HashMap<_, _, BuildHasherDefault<XxHash64>> = Default::default();
    for result in rdr.records() {
        let record = result?;
        map.insert(record[key_side].to_owned(), record[1-key_side].to_owned());
        // eprintln!("Loaded {}", record[key_side].to_owned());
    }
    Ok(map)
}