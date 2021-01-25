use std::path::PathBuf;
use structopt::StructOpt;

pub mod releases;
pub mod normalize;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Subcommand {
    /// Find the max release for each urs entry.
    MaxRelease {
        #[structopt(parse(from_os_str))]
        /// Filename of the raw json file, '-' means stdin.
        filename: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    /// Take the output of the kv lookup and turn it into JSON suitable for the pipeline
    Normalize {
        #[structopt(parse(from_os_str))]
        /// Filename of the raw json file, '-' means stdin.
        filename: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    /// Select all xrefs to use. Xrefs will be selected if they have a newer release than the
    /// previous one, if they have never been precomputed or if they are new.
    Select {
        #[structopt(parse(from_os_str))]
        /// Filename of the xref entries, must be sorted
        xrefs: PathBuf,
        
        #[structopt(parse(from_os_str))]
        /// Filename of the known entries, must be sorted
        known: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    }
}

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
struct Opt {
    /// Set the logging option, more is more verbose.
    #[structopt(short = "v", long = "verbose", parse(from_occurrences))]
    verbose: u32,

    #[structopt(subcommand)]
    command: Subcommand,
}

fn main() -> anyhow::Result<()> {
    let opt = Opt::from_args();

    let level = match opt.verbose {
        0 => simplelog::LevelFilter::Warn,
        1 => simplelog::LevelFilter::Info,
        2 => simplelog::LevelFilter::Debug,
        _ => simplelog::LevelFilter::Trace,
    };
    simplelog::TermLogger::init(
        level,
        simplelog::Config::default(),
        simplelog::TerminalMode::Stderr,
    )
    .unwrap_or_else(|_| eprintln!("Failed to create logger, ignore"));

    match opt.command {
        Subcommand::MaxRelease { filename, output } => {
            releases::write_max(&filename, &output)?;
        }
        Subcommand::Normalize { filename, output } => {
            normalize::write(&filename, &output)?;
        }
        Subcommand::Select { xrefs, known, output } => {
            releases::select_new(&xrefs, &known, &output)?;
        }
    }

    Ok(())
}
