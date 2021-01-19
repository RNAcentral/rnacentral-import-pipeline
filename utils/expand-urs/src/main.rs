use std::path::{
    Path,
    PathBuf,
};

use structopt::StructOpt;

use anyhow::{
    Context,
    Result,
};

use rnc_core::{
    containers::urs_taxid::UrsTaxidMapping,
    urs::Urs,
};

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Subcommand {
    /// This is a tool to process a file of JSON objects and expand their URS entry to
    /// urs_taxid entries. Each object must contain a 'urs' field which contains the URS
    /// to expand. This will then produce an object with an 'id' field for each
    /// urs_taxid with the same urs.
    Json {
        /// Name of the field in the JSON object to expand.
        #[structopt(long, short, default_value = "urs")]
        field_name: String,

        /// A file where each line is a urs_taxid, which are all urs_taxids that need to
        /// be output. Duplicate will be treated as single entry.
        #[structopt(parse(from_os_str))]
        active_file: PathBuf,

        /// A file ('-' means stdin) where each line is a valid json object, which
        /// contains a key 'urs' that is a URS that exists in RNAcentral. If the
        /// UPI is not in the active_file the object will not be written and a
        /// warning logged.
        #[structopt(parse(from_os_str))]
        filename: PathBuf,

        /// File to output to, '-' means stdout.
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },

    Text {
        /// A file where each line is a urs_taxid, which are all urs_taxids that need to
        /// be output. Duplicate will be treated as single entry.
        #[structopt(parse(from_os_str))]
        active_file: PathBuf,

        /// A file ('-' means stdin) where each line is a valid json object, which
        /// contains a key 'urs' that is a URS that exists in RNAcentral. If the
        /// UPI is not in the active_file the object will not be written and a
        /// warning logged.
        #[structopt(parse(from_os_str))]
        filename: PathBuf,

        /// File to output to, '-' means stdout.
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },
}

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
/// This is a tool to expand URS to all possible urs_taxid pairs. This works for a variety
/// of formats.
struct Opt {
    #[structopt(subcommand)]
    command: Subcommand,
}

fn expand_json(active_file: &Path, filename: &Path, field_name: &str, output: &Path) -> Result<()> {
    let mut input = rnc_utils::buf_reader(filename)?;
    let mut output = rnc_utils::buf_writer(output)?;

    let container = UrsTaxidMapping::from_urs_file(active_file)?;
    let mut buf = String::new();
    loop {
        match input.read_line(&mut buf)? {
            0 => break,
            _ => {
                let line = buf.replace("\\\\", "\\");
                let mut json: serde_json::Value = serde_json::from_str(&line)
                    .with_context(|| format!("Cannot parse JSON object {}", &line))?;

                if let Some(m) = json.as_object_mut() {
                    if let Some(serde_json::Value::String(raw_urs)) = m.get(field_name) {
                        let urs: Urs = raw_urs
                            .parse()
                            .with_context(|| format!("Failed to parse URS id {}", &raw_urs))?;
                        let urs_taxids = container.urs_taxids(&urs);
                        if urs_taxids.len() == 0 {
                            log::warn!("No active URS_taxid found for: {}", &raw_urs);
                        }

                        for urs_taxid in urs_taxids {
                            m.insert(
                                String::from("id"),
                                serde_json::Value::String(urs_taxid.to_string()),
                            );
                            serde_json::to_writer(&mut output, &m)?;
                            writeln!(&mut output)?;
                        }
                    }
                }

                buf.clear();
            },
        }
    }

    Ok(())
}

fn expand_text(active_file: &Path, filename: &Path, output: &Path) -> Result<()> {
    let mut input = rnc_utils::buf_reader(filename)?;
    let mut output = rnc_utils::buf_writer(output)?;

    let container = UrsTaxidMapping::from_urs_file(active_file)?;
    let mut buf = String::new();
    loop {
        match input.read_line(&mut buf)? {
            0 => break,
            _ => {
                let urs: Urs =
                    buf.parse().with_context(|| format!("Failed to parse URS id {}", &buf))?;
                let urs_taxids = container.urs_taxids(&urs);
                if urs_taxids.len() == 0 {
                    log::warn!("No active URS_taxid found for: {}", &buf);
                }

                for urs_taxid in urs_taxids {
                    writeln!(&mut output, "{}", urs_taxid.to_string())?;
                }

                buf.clear();
            },
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let opt = Opt::from_args();
    match opt.command {
        Subcommand::Json {
            field_name,
            active_file,
            filename,
            output,
        } => expand_json(&active_file, &filename, &field_name, &output)?,
        Subcommand::Text {
            active_file,
            filename,
            output,
        } => expand_text(&active_file, &filename, &output)?,
    }
    Ok(())
}
