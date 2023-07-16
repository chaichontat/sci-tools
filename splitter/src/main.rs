#![feature(portable_simd)]
#![feature(ascii_char)]
pub mod simdstuffs;
use clap::Parser;

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::simd::{u8x32, SimdPartialEq};

use crate::simdstuffs::mismatch_count;
use csv::ReaderBuilder;
use lazy_static::lazy_static;
use needletail::parse_fastx_file;

const MAX_SUBSEQ_LEN: usize = 32;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    fastq_file: String,
    tag_file: String,
    #[arg(short, long)]
    output: String,
    #[arg(short, long)]
    index: usize,
    #[arg(short, long)]
    tol: Option<i8>,
}

fn main() {
    let args = Args::parse();

    let fastq_file = &args.fastq_file;
    let tag_file = &args.tag_file;
    let output_file = &args.output;
    let tol = args.tol.unwrap_or(0);
    if !(0..=32).contains(&tol) {
        panic!("Tolerance must be between 0 and 32");
    }

    let mut tags = Vec::new();
    let mut r = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(tag_file)
        .unwrap();

    let start = args.index;
    let mut end = 0;

    for (i, result) in r.records().enumerate() {
        let record = result.unwrap();
        if i == 0 {
            end = start + record.get(0).unwrap().len();
        } else if record.get(0).unwrap().len() + start != end {
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

    let mut output = BufWriter::new(File::create(output_file).unwrap());
    let mut reader = parse_fastx_file(fastq_file).expect("valid path/file");
    let mut count = 0;
    let mut unk = 0;

    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        // println!("{:?}", str::from_utf8(seqrec.raw_seq()[start..end].as_ascii().unwrap().as_bytes()));
        let sub_seq = encode_sequence(seqrec.raw_seq()[start..end].as_ascii().unwrap().as_bytes());
        'found: {
            for (sequence, name) in &tags {
                if mismatch_count(&sub_seq, sequence, tol) {
                    writeln!(output, "{}", name).unwrap();
                    break 'found;
                }
            }
            unk += 1;
            writeln!(output, "unk").unwrap();
        }

        count += 1;
        if count % 1000000 == 0 {
            println!(
                "Processed {} records. Unknown {:.4}.",
                count,
                unk as f32 / count as f32
            );
            if unk as f32 / count as f32 > 1. {
                println!("Too many unknowns");
                break;
            }
        }
        if count > 100000 {
            break;
        }
    }
    output.flush().unwrap();
}

fn encode_sequence(sequence: &[u8]) -> [u8; MAX_SUBSEQ_LEN] {
    let mut encoded = [0; MAX_SUBSEQ_LEN];
    for (i, &nucleotide) in sequence.iter().enumerate() {
        encoded[i] = nucleotide
    }
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
