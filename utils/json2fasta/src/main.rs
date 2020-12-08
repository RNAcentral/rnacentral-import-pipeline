use std::error::Error;
use std::path::PathBuf;

extern crate log;

use bio::io::fasta;

use structopt::StructOpt;

use rnc_core::json_sequence::Sequence;
use rnc_core::nhmmer::valid_sequence;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
struct Opt {
    /// This will select only sequence which have known easel/infernal characters. This limits things
    /// to matching ACGUN.
    #[structopt(short, long)]
    only_valid_easel: bool,

    /// The name of the file to read from, using '-' means stdin.
    #[structopt(parse(from_os_str))]
    raw: PathBuf,

    /// The name of the file to write to, using - means stdout.
    #[structopt(parse(from_os_str))]
    output: PathBuf,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();
    let mut reader = rnc_utils::buf_reader(&opt.raw)?;
    let output = rnc_utils::buf_writer(&opt.output)?;
    let mut writer = fasta::Writer::new(output);
    let mut buf = String::new();

    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let cleaned = buf.replace("\\\\", "\\");
                let data: Sequence = serde_json::from_str(&cleaned)?;
                if !opt.only_valid_easel || valid_sequence(&data.sequence) {
                    writer.write_record(&data.into())?;
                }
                buf.clear();
            }
        }
    }

    Ok(())
}
