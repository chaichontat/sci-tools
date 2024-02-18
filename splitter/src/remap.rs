use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;
use clap::Parser;
use scitools::extract_idx;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = "Simple remapper")]
struct Args {
    /// Path to tag file (TSV format)
    tag_file: String,
}



fn main() -> io::Result<()> {
    let args = Args::parse();
    let tag_file = &args.tag_file;
    let hash = extract_idx(Path::new(tag_file), 0).expect("Hash failed");

    eprintln!("Loaded {} records", hash.len());
    let mut buf = String::new();
    let mut stdin = io::stdin().lock();
    // let mut output = BufWriter::new(std::io::stdout());

    let mut i = 0;
    let unk = &"0".to_string();

    while let Ok(line) = stdin.read_line(&mut buf) {
        if line == 0 {
            break;
        }
        let curr = buf.split("\t").next().expect("No split");
        let record = hash.get(curr).unwrap_or(unk);
        println!("{}", record);
        i += 1;
        if i % 1000000 == 0 {
            eprintln!("Processed {} records", i);
        }
        buf.clear();
    }


    Ok(())
}