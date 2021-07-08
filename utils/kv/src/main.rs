use std::{
    error::Error,
    path::PathBuf,
};
extern crate log;
use structopt::StructOpt;

pub mod store;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Subcommand {
    /// Index a JSON document.
    Index {
        #[structopt(short, long, default_value = "1000000")]
        commit_size: usize,

        /// Type of data being indexed, eg, secondary_structure, hits, etc
        data_type: String,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw json file, '-' means stdin.
        filename: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to store the index in.
        output: PathBuf,
    },

    /// Index a file of JSON objects where the objects are sorted by id.
    SortedIndex {
        #[structopt(short, long, default_value = "1000000")]
        commit_size: usize,

        /// Type of data being indexed, eg, secondary_structure, hits, etc
        data_type: String,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw json file, '-' means stdin.
        filename: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to store the index in.
        output: PathBuf,
    },

    /// Given a file where each line is an id to extract, extract all values for it.
    Lookup {
        /// This will cause the program to emit a warning instead of fail if a requested
        /// key has no data.
        #[structopt(short, long)]
        allow_missing: bool,

        /// Filename of the database file
        #[structopt(parse(from_os_str))]
        cache: PathBuf,

        /// Filename of the id files, '-' means stdin.
        #[structopt(parse(from_os_str))]
        filename: PathBuf,

        /// Filename to write the data to, '-' means stdout.
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },
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

fn main() -> Result<(), Box<dyn Error>> {
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
        Subcommand::Index {
            commit_size,
            data_type,
            filename,
            output,
        } => {
            let mut spec = store::Spec::new(&output);
            spec.set_commit_size(commit_size);
            store::index(&spec, &data_type, &filename)?
        },
        Subcommand::SortedIndex {
            commit_size,
            data_type,
            filename,
            output,
        } => {
            let mut spec = store::Spec::new(&output);
            spec.set_commit_size(commit_size);
            store::sorted_index(&spec, &data_type, &filename)?
        },
        Subcommand::Lookup {
            allow_missing,
            cache,
            filename,
            output,
        } => {
            let mut spec = store::Spec::new(&cache);
            spec.set_allow_missing(allow_missing);
            store::lookup(&spec, &filename, &output)?
        },
    };

    Ok(())
}
