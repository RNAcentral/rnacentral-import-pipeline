use clap::Parser;
use std::fs;
use polars::datatypes::DataType::Int64;
use polars::prelude::*;

#[derive(Parser, Debug)]
#[clap(author = "Andrew Green", version, about)]
struct Args {
    input: String,
    output: String,
}

/// Take the standard RNAcentral bed file and expand it out so it has one line per exon
fn main() -> Result<(), PolarsError> {
    let cli = Args::parse();

    let original_bed =
        LazyCsvReader::new(cli.input).has_header(false).with_delimiter(b'\t').finish().unwrap();

    // Stay lazy as long as possible to minimise memory use
    let original_bed = original_bed
        .with_columns([col("column_11").str().split(","), col("column_12").str().split(",")])
        .explode([col("column_11"), col("column_12")])
        .with_columns([col("column_11").cast(Int64), col("column_12").cast(Int64)])
        .with_columns([(col("column_2") + col("column_12")).alias("column_2"), (col("column_2") + col("column_12") + col("column_11")).alias("column_3")] )
        .select([col("column_1"), col("column_2"), col("column_3"), col("column_4"), col("column_5"), col("column_6")])
        .sort_by_exprs(vec![col("column_1"), col("column_2")], vec![false,false], false);

    let mut output_file = fs::File::create(&cli.output)?;
    let mut writer = CsvWriter::new(&mut output_file).has_header(false).with_delimiter(b'\t');

    let mut expanded_bed = original_bed.collect()?;

    writer.finish(&mut expanded_bed)

}
