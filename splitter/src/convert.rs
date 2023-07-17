#![feature(portable_simd)]
#![feature(ascii_char)]

use clap::Parser;

use needletail::{parse_fastx_file, parse_fastx_stdin};
use parquet::{
    basic::Compression::ZSTD, basic::ZstdLevel, data_type::ByteArrayType,
    file::properties::WriterProperties, schema::parser::parse_message_type,
};
use std::{fs::File, path::Path};

use parquet::{data_type::ByteArray, file::writer::SerializedFileWriter};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    fastq_file: String,
    output: String,
}

const ROW_GROUP_SIZE: usize = (2 as usize).pow(25);
use std::{fs, sync::Arc};

fn create_file(p: &str) -> SerializedFileWriter<File> {
    let file = fs::File::create(&Path::new(p)).unwrap();
    let message_type = "
        message schema {
            REQUIRED BINARY name(UTF8);
            REQUIRED BINARY sequence(UTF8);
            REQUIRED BINARY quality(UTF8);
        }
    ";
    let schema = Arc::new(parse_message_type(message_type).unwrap());
    let props = Arc::new(
        WriterProperties::builder()
            .set_compression(ZSTD(ZstdLevel::try_new(4).unwrap()))
            .build(),
    );
    SerializedFileWriter::new(file, schema, props).unwrap()
}

fn main() {
    let args = Args::parse();

    let mut reader = match (&args.fastq_file).as_str() {
        "-" => parse_fastx_stdin().expect("valid stdin"),
        _ => parse_fastx_file(args.fastq_file).expect("valid path/file"),
    };

    let mut i = 0;
    let mut total = 0;
    let mut writer = create_file(&args.output);

    let mut temp: Vec<Vec<Vec<u8>>> = vec![Vec::new(), Vec::new(), Vec::new()];

    while let Some(record) = reader.next() {
        if i & ROW_GROUP_SIZE != 0 {
            write_row_group(&mut writer, &temp);
            i = 0;
            for element in &mut temp {
                element.clear();
            }
        }

        let seqrec = record.expect("invalid record");
        temp[0].push(seqrec.id().to_vec());
        temp[1].push(seqrec.seq().to_vec());
        temp[2].push(seqrec.qual().unwrap().to_vec());
        i += 1;
        total += 1;

        if total % (2 << 20) == 0 {
            println!("{} records processed", total);
        }
    }

    if i > 0 {
        write_row_group(&mut writer, &temp);
    }
    writer.close().unwrap();
}

fn write_row_group(writer: &mut SerializedFileWriter<File>, data: &[Vec<Vec<u8>>]) {
    assert!(data.len() > 0);
    let mut row_group_writer = writer.next_row_group().unwrap();
    let m = data[0].len();

    for i in 0..m {
        if let Some(mut col_writer) = row_group_writer.next_column().unwrap() {
            // let values = ByteArray::from(value);
            col_writer
                .typed::<ByteArrayType>()
                .write_batch(
                    &data[i]
                        .iter()
                        .map(|x| ByteArray::from(x.as_slice()))
                        .collect::<Vec<parquet::data_type::ByteArray>>(),
                    None,
                    None,
                )
                .unwrap();
            col_writer.close().unwrap();
        }
    }
    row_group_writer.close().unwrap();
}

macro_rules! vec_of_strings {
    ($($x:expr),*) => (vec![$($x.to_string()),*]);
}

#[test]
fn sample_test() {
    use parquet::file::properties::WriterProperties;
    use parquet::schema::parser::parse_message_type;
    use std::path::Path;
    use std::{fs, sync::Arc};

    let path = Path::new("./sample.parquet");
    let message_type = "
        message schema {
            REQUIRED BINARY name(UTF8);
            REQUIRED BINARY sequence(UTF8);
            REQUIRED BINARY quality(UTF8);
        }
    ";
    let schema = Arc::new(parse_message_type(message_type).unwrap());
    let props = Arc::new(WriterProperties::builder().build());
    let file = fs::File::create(&path).unwrap();

    let dummy: Vec<u8> = vec![0, 0, 0];
    let data = vec![vec![dummy.clone(); 10]; 3];

    let mut writer = SerializedFileWriter::new(file, schema, props).unwrap();
    write_row_group(&mut writer, &data);
    writer.close().unwrap();
}
