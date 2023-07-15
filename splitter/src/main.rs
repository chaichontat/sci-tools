#![feature(portable_simd)]
#![feature(ascii_char)]
use std::env;
use std::str;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::simd::{u8x32, SimdPartialEq};

use csv::ReaderBuilder;
use needletail::parse_fastx_file;

const MAX_SUBSEQ_LEN: usize = 32;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!(
            "Usage: {} <fastq_file> <location> <tag_file> <output_file>",
            args[0]
        );
        std::process::exit(1);
    }

    let fastq_file = &args[1];
    let location_str = &args[2];
    let tag_file = &args[3];
    let output_file = &args[4];
    let tol =  {
        if args.len() < 5 + 1 {
            2
         } else {
            args[5].parse::<i8>().expect("tolerance must be number")
        }
    } as i8;

    // let mut loc_iter = location_str.split(',');

    // let end: usize = loc_iter.next().unwrap().parse().unwrap();

    let mut tags = Vec::new();
    let mut r = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(tag_file)
        .unwrap();

    let start: usize = location_str.parse().expect("start must be number");
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
                if is_one_mismatch_or_less(&sub_seq, sequence, tol) {
                    writeln!(output, "{}", name).unwrap();
                    break 'found;
                }
            }
            unk += 1;
            writeln!(output, "unk").unwrap();
        }

        count += 1;
        if count % 1000000 == 0 {
            println!("Processed {} records. Unknown {:.4}.", count, unk as f32 / count as f32);
            if unk as f32 / count as f32  > 1. {
                println!("Too many unknowns");
                break;
            }
        }
        // if count > 33 {
        //     break;
        // }
    }
}

fn encode_sequence(sequence: &[u8]) -> [u8; MAX_SUBSEQ_LEN] {
    let mut encoded = [0; MAX_SUBSEQ_LEN];
    for (i, &nucleotide) in sequence.iter().enumerate() {
        encoded[i] |= encode_nucleotide(nucleotide);
    }
    encoded
}

fn encode_nucleotide(nucleotide: u8) -> u8 {
    let x = match nucleotide.to_ascii_uppercase() {
        b'A' => 0b00,
        b'T' => 0b01,
        b'C' => 0b10,
        b'G' => 0b11,
        b'N' => 0b100,
        x => x,
    };
    if x == nucleotide {
        println!("Unknown nucleotide: {}", nucleotide);
    }
    x
}

fn is_one_mismatch_or_less(
    sub_seq: &[u8; MAX_SUBSEQ_LEN],
    sequence: &[u8; MAX_SUBSEQ_LEN],
    tol: i8,
) -> bool {
    let simd_sub_seq = u8x32::from_slice(sub_seq);
    let simd_sequence = u8x32::from_slice(sequence);
    let simd_result = simd_sub_seq.simd_eq(simd_sequence);
    let mismatch_count = simd_result.to_int().to_array();
    let res = mismatch_count.iter().sum::<i8>();
    res + 32 <= tol
}

// fn main() {
//     let a = u8x32::splat(10.0);
//     let b = u8x32::from_array([1.0, 2.0, 3.0, 4.0]);
//     println!("{:?}", a + b);
// }
