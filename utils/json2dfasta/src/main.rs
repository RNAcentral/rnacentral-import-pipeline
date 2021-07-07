use std::path::PathBuf;

use anyhow::{
    Context,
    Result,
};

use serde::Deserialize;

use structopt::StructOpt;

#[derive(Debug, Deserialize)]
struct Entry {
    id: String,
    sequence: String,

    #[serde(rename = "secondary-structure")]
    secondary_structure: String,
}

/// This is a utility script to convert JSON objects to FASTA files with dot-bracket base
/// pairs. The JSON files should be formatted like:
///
/// { "id": $urs, "sequence": $sequence, "secondary-structure": $dot-bracket }
///
/// and will be formatted to produce files like:
///
/// >$urs \n
/// $sequence \n
/// $dot-bracket \n
///
/// That is to say a FASTA file with the ID as the URS, the sequence on one line, and the
/// base pairs on another.
///
/// This is a useful format for people that want to download our entire secondary
/// structure data.
#[derive(Debug, StructOpt)]
struct Cli {
    #[structopt(parse(from_os_str))]
    /// The input file name, may be '-' for stdin.
    input: PathBuf,

    #[structopt(parse(from_os_str))]
    /// The output file name, may be '-' for stdout.
    output: PathBuf,
}

fn main() -> Result<()> {
    let opt = Cli::from_args();

    let mut reader = rnc_utils::buf_reader(&opt.input)?;
    let mut writer = rnc_utils::buf_writer(&opt.output)?;

    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let entry: Entry = serde_json::from_str(&buf)
                    .with_context(|| format!("Could not parse line: {}", &buf))?;
                writeln!(&mut writer, ">{}", entry.id).with_context(|| {
                    format!("Failed writing id from {:?} to {:?}", &entry, &opt.output)
                })?;
                writeln!(&mut writer, "{}", entry.sequence).with_context(|| {
                    format!("Failed writing sequence from {:?} to {:?}", &entry, &opt.output)
                })?;
                writeln!(&mut writer, "{}", entry.secondary_structure).with_context(|| {
                    format!("Failed writing pairs from {:?} to {:?}", &entry, &opt.output)
                })?;
                buf.clear();
            },
        }
    }

    Ok(())
}
