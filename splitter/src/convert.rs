#![feature(generators, generator_trait)]
use clap::{Parser, Subcommand};
use io::to_parquet;

use crate::io::ParquetReader;

pub mod io;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: SubCommand,
}

#[derive(Debug, Subcommand)]
enum SubCommand {
    Convert {
        fastq_file: String,
        output: String,
    },
    Read {
        parquet_file: String,
        output: String,
    },
}

fn main() {
    let args = Args::parse();

    match &args.command {
        SubCommand::Convert { fastq_file, output } => {
            to_parquet(fastq_file, output);
        }
        SubCommand::Read {
            parquet_file,
            output,
        } => {
            // read_parquet(parquet_file);
        }
    }
}

#[test]
fn sample_test() {
    use parquet::file::properties::WriterProperties;
    use parquet::schema::parser::parse_message_type;
    use std::path::Path;
    use std::{fs, sync::Arc};
    let mut i = 0;

    for c in ParquetReader::new("out.parquet") {
        i += 1;
        if i % 1000000 == 0 {
            println!("{} records processed", i);
        }
    }
}
