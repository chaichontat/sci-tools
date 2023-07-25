use needletail::{parse_fastx_file, parse_fastx_stdin};
use parquet::{
    basic::Compression::ZSTD,
    basic::ZstdLevel,
    basic::{ConvertedType, Repetition, Type as PhysicalType},
    data_type::{ByteArrayType, FixedLenByteArray, FixedLenByteArrayType},
    file::{
        properties::WriterProperties,
        reader::{FileReader, RowGroupReader, SerializedFileReader},
    },
    record::Field,
    schema::{parser::parse_message_type, types::Type},
};
use parquet::{data_type::ByteArray, file::writer::SerializedFileWriter};
use std::{fs, sync::Arc};
use std::{fs::File, path::Path};

const ROW_GROUP_SIZE: usize = (2 as usize).pow(23);

fn create_file(p: &str, length: usize) -> SerializedFileWriter<File> {
    let file = fs::File::create(&Path::new(p)).unwrap();
    let message_type = format!(
        "
        message schema {{
            REQUIRED BINARY name(UTF8);
            REQUIRED FIXED_LEN_BYTE_ARRAY({}) sequence;
            REQUIRED FIXED_LEN_BYTE_ARRAY({}) quality;
        }}
    ",
        length, length
    );

    let schema = Arc::new(parse_message_type(&message_type).unwrap());
    let props = Arc::new(
        WriterProperties::builder()
            .set_compression(ZSTD(ZstdLevel::try_new(4).unwrap()))
            .build(),
    );
    SerializedFileWriter::new(file, schema, props).unwrap()
}

fn write_row_group(writer: &mut SerializedFileWriter<File>, data: &mut [Vec<Vec<u8>>]) {
    assert!(data.len() > 0);
    let mut row_group_writer = writer.next_row_group().unwrap();
    let m = data[0].len();

    for i in 0..m {
        if let Some(mut col_writer) = row_group_writer.next_column().unwrap() {
            if i == 0 {
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
            } else {
                col_writer
                    .typed::<FixedLenByteArrayType>()
                    .write_batch(
                        &data[i]
                            .iter_mut()
                            .map(|x| FixedLenByteArray::from(std::mem::replace(x, Vec::new())))
                            .collect::<Vec<parquet::data_type::FixedLenByteArray>>(),
                        None,
                        None,
                    )
                    .unwrap();
            }
            col_writer.close().unwrap();
        }
    }
    row_group_writer.close().unwrap();
}

pub struct ParquetReader {
    reader: SerializedFileReader<File>,
    curr_idx: usize,
    projection: Type,
    curr_rows: Vec<Vec<u8>>,
    next_row_group: usize,
    n_row_group: usize,
}

impl ParquetReader {
    pub fn new(path: &str) -> Self {
        let file = File::open(&Path::new(path)).expect("Couldn't open parquet file");
        let reader = SerializedFileReader::new(file).unwrap();
        let metadata = reader.metadata();
        let length = metadata.row_group(0).columns()[1]
            .column_descr()
            .type_length();

        let field_a = Type::primitive_type_builder("sequence", PhysicalType::FIXED_LEN_BYTE_ARRAY)
            .with_length(length)
            .with_repetition(Repetition::REQUIRED)
            .build()
            .unwrap();

        let projection = Type::group_type_builder("schema")
            .with_fields(&mut vec![Arc::new(field_a)])
            .build()
            .unwrap();

        let n = reader.num_row_groups();

        ParquetReader {
            reader,
            projection,
            curr_idx: 0,
            curr_rows: Vec::new(),
            next_row_group: 0,
            n_row_group: n,
        }
    }

    fn _get(&self) -> Option<Vec<u8>> {
        let out = Some(self.curr_rows[self.curr_idx].clone());
        return out;
    }
}

impl Iterator for ParquetReader {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr_rows.len() > 0 && self.curr_idx < self.curr_rows.len() - 1 {
            let out = self._get();
            self.curr_idx += 1;
            return out;
        }
        if self.curr_rows.len() == 0
            || (self.curr_idx == self.curr_rows.len() - 1 && self.next_row_group < self.n_row_group)
        {
            // Go to next row group
            let rg_reader = self.reader.get_row_group(self.next_row_group).unwrap();
            self.curr_rows = row_group_reader(&*rg_reader, self.projection.clone());
            self.curr_idx = 0;
            self.next_row_group += 1;

            let out = self._get();
            self.curr_idx += 1;
            return out;
        }
        None
    }
}

// fn read_parquet(path: &str, function: fn(&[Vec<u8>])) {
//     let file = File::open(&Path::new(path)).expect("Couldn't open parquet file");
//     let reader = SerializedFileReader::new(file).unwrap();
//     let metadata = reader.metadata();
//     let length = metadata.row_group(0).columns()[1]
//         .column_descr()
//         .type_length();

//     let field_a = Type::primitive_type_builder("sequence", PhysicalType::FIXED_LEN_BYTE_ARRAY)
//         .with_length(length)
//         .with_repetition(Repetition::REQUIRED)
//         .build()
//         .unwrap();

//     let projection = Type::group_type_builder("schema")
//         .with_fields(&mut vec![Arc::new(field_a)])
//         .build()
//         .unwrap();

//     let n = reader.num_row_groups();
//     for i in 0..n {
//         let rg_reader = reader.get_row_group(i).unwrap();
//         let data = row_group_reader(&*rg_reader, projection.clone());
//         function(&data[..])
//     }
// }

fn row_group_reader<'a>(reader: &'a dyn RowGroupReader, projection: Type) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    if let Ok(it) = reader.get_row_iter(Some(projection)) {
        for row in it {
            let r = row.unwrap();
            let field = r.get_column_iter().next().unwrap().1;
            if let Field::Bytes(x) = field {
                out.push(x.as_ref().to_vec());
                // let y: &[u8] = x.as_ref();
            }
        }
    }
    out
}

pub fn to_parquet(fastq_file: &str, output: &str) {
    let mut reader = match fastq_file {
        "-" => parse_fastx_stdin().expect("valid stdin"),
        _ => parse_fastx_file(fastq_file).expect("valid path/file"),
    };

    let mut i = 0;
    let mut total = 0;
    let mut writer: Option<SerializedFileWriter<File>> = None;

    let mut temp: Vec<Vec<Vec<u8>>> = vec![Vec::new(), Vec::new(), Vec::new()];

    while let Some(record) = reader.next() {
        if total == 0 {
            writer = Some(create_file(output, record.as_ref().unwrap().seq().len()));
        }
        if i & ROW_GROUP_SIZE != 0 {
            write_row_group(&mut writer.as_mut().unwrap(), &mut temp);
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

    if i > 0 && writer.is_some() {
        write_row_group(&mut writer.as_mut().unwrap(), &mut temp);
    }

    if writer.is_some() {
        writer.unwrap().close().unwrap();
    }
}
