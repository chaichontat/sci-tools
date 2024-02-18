#![feature(portable_simd)]
#![feature(ascii_char)]

pub mod io;
pub mod simdstuffs;
use std::path::Path;
use clap::Parser;
use scitools::FastqSequence;

use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::str;

use crate::simdstuffs::mismatch_count;
use csv::ReaderBuilder;
use needletail::{parse_fastx_file, parse_fastx_stdin, FastxReader};

const MAX_SUBSEQ_LEN: usize = 32;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = "Read demuxer for sci-rna-seq3.")]
struct Args {
    /// Path to read 1 FASTQ file
    fastq_file: String,
    /// Path to tag file (TSV format)
    tag_file: String,
    /// Index in read 1 to start matching (default=0)
    #[arg(short, long)]
    index: usize,
    /// Mismatch tolerance (default=1)
    #[arg(short, long)]
    tol: Option<i8>,
    #[arg(short, long)]
    /// Whether to read for shifts
    shift_file: Option<String>,
}

fn parse_tag(tag_file: &str) -> (Vec<([u8; MAX_SUBSEQ_LEN], String)>, usize) {
    let mut tags = Vec::new();
    let mut r = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(tag_file)
        .unwrap();

    let mut len = 0;

    for (i, result) in r.records().enumerate() {
        let record = result.unwrap();
        if i == 0 {
            len = record.get(0).unwrap().len();
        } else if record.get(0).unwrap().len() != len {
            panic!("All sequences must be same length");
        }

        let sequence = encode_sequence(
            record
                .get(0)
                .unwrap()
                .to_ascii_uppercase()
                .as_ascii()
                .unwrap()
                .as_bytes(),
        );
        let name = record.get(1).unwrap().to_string();
        tags.push((sequence, name));
    }
    (tags, len)
}



fn main() {
    let args = Args::parse();

    let fastq_file = &args.fastq_file;
    let tag_file = &args.tag_file;
    let shift_file = &args.shift_file;

    let mut output = BufWriter::new(std::io::stdout());

    let mut shift_reader: Box<dyn Iterator<Item = usize>> = if let Some(shift_file) = shift_file {
        eprintln!("Reading shifts");
        let path = Path::new(shift_file);
        let file = File::open(path).unwrap();
        Box::new(std::io::BufReader::new(file).lines().map(|x| x.unwrap().parse::<usize>().unwrap()))
    } else {
        // Iterator that outputs a constant value
        Box::new(std::iter::repeat(args.index))
    };


    let tol = args.tol.unwrap_or(0);
    if !(0..=32).contains(&tol) {
        panic!("Tolerance must be between 0 and 32");
    }
    let (tags, len) = parse_tag(tag_file);
    let start = args.index;
    let end = start + len;

    // let mut output = BufWriter::new(File::create(output_file).unwrap());
    // let mut output be the stdout


    let mut reader: Box<dyn Iterator<Item = (Vec<u8>, Vec<u8>)>> = match fastq_file.as_str() {
        "-" => Box::new(FastqSequence {
            reader: parse_fastx_stdin().expect("valid stdin"),
        }),
        x if x.ends_with(".fastq") => Box::new(FastqSequence {
            reader: parse_fastx_file(fastq_file).expect("valid path/file"),
        }),
        _ => panic!("Invalid file type")
    };
    let mut count = 0;
    let mut unk = 0;

    while let Some((seq, qual)) = reader.next() {
        // println!("{:?}", str::from_utf8(seqrec.raw_seq()[start..end].as_ascii().unwrap().as_bytes()));
        let s = shift_reader.next().unwrap();
        let sub_seq = encode_sequence(&seq[start+s..end+s]);
        'found: {
            let qual = str::from_utf8(&qual[start+s..end+s]).unwrap();
            for (sequence, name) in &tags {
                if mismatch_count(&sub_seq, sequence, tol) {
                    writeln!(output, "{}\t{}", name, qual).unwrap();
                    break 'found;
                }
            }
            unk += 1;
            writeln!(output, "unk\t{}", qual).unwrap();
        }

        count += 1;
        if count % 1000000 == 0 {
            eprintln!(
                "Processed {} records. Unknown {:.4}.",
                count,
                unk as f32 / count as f32
            );
            if unk as f32 / count as f32 > 1. {
                eprintln!("Too many unknowns");
                break;
            }
        }
    }
    output.flush().unwrap();
}

fn encode_sequence(sequence: &[u8]) -> [u8; MAX_SUBSEQ_LEN] {
    let mut encoded = [0; MAX_SUBSEQ_LEN];
    encoded[0..sequence.len()].clone_from_slice(sequence);
    encoded
}

// fn encode_nucleotide(nucleotide: u8) -> u8 {
//     let x = match nucleotide.to_ascii_uppercase() {
//         b'A' => 0b00,
//         b'T' => 0b01,
//         b'C' => 0b10,
//         b'G' => 0b11,
//         b'N' => 0b100,
//         x => x,
//     };
//     if x == nucleotide {
//         println!("Unknown nucleotide: {}", nucleotide);
//     }
//     x
// }

// lazy_static! {
// static ref MASK: Simd<u8, 32> = u8x32::splat(32);
// };

// fn main() {
//     let a = u8x32::splat(10.0);
//     let b = u8x32::from_array([1.0, 2.0, 3.0, 4.0]);
//     println!("{:?}", a + b);
// }
