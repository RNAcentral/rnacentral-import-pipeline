use std::{
    path::PathBuf,
    str::FromStr,
};
use structopt::StructOpt;

pub mod accessions;
pub mod metadata;
pub mod normalize;
pub mod releases;

#[derive(Debug)]
enum Groupable {
    Basic,
    Coordinates,
    Previous,
    R2dtHits,
    RfamHits,
}

impl FromStr for Groupable {
    type Err = String;

    fn from_str(raw: &str) -> Result<Self, Self::Err> {
        match raw {
            "basic" => Ok(Self::Basic),
            "coordinates" => Ok(Self::Coordinates),
            "previous" => Ok(Self::Previous),
            "r2dt-hits" => Ok(Self::R2dtHits),
            "r2dt_hits" => Ok(Self::R2dtHits),
            "rfam-hits" => Ok(Self::RfamHits),
            "rfam_hits" => Ok(Self::RfamHits),
            unknown => Err(format!("Unknown name {}", unknown)),
        }
    }
}

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum MetadataCommand {
    /// Merge each distinct type of metadata into single merged entries. THis assumes that
    /// all files are sorted in the same way. Additionally, the xref file must have at
    /// least one entry for all possible urs_taxids.
    Merge {
        #[structopt(parse(from_os_str))]
        /// Filename of the raw xref file
        basic: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw coordinates file
        coordinates: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw rfam hits file
        rfam_hits: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw r2dt hits file
        r2dt_hits: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the raw previous file
        previous: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    Group {
        /// Type of data to group
        #[structopt(case_insensitive = true)]
        data_type: Groupable,

        /// Filename to read the results from, '-' means stdin
        #[structopt(parse(from_os_str))]
        path: PathBuf,

        /// The maximum count of the entries.
        max_count: usize,

        /// Filename to write the results to, '-' means stdout
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },
}

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
        /// Filename of the raw accession json file, '-' means stdin.
        accessions: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the merged metadata file
        metadata: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    /// Commands dealing with metadata.
    Metadata {
        #[structopt(subcommand)]
        command: MetadataCommand,
    },

    GroupAccessions {
        /// Filename to read the results from, '-' means stdin
        #[structopt(parse(from_os_str))]
        path: PathBuf,

        /// The min index of these accessions.
        min_index: usize,

        /// The max index of these accessions.
        max_index: usize,

        /// Filename to write the results to, '-' means stdout
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },

    /// Select all xrefs to use. Xrefs will be selected if they have a newer release than
    /// the previous one, if they have never been precomputed or if they are new.
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
        Subcommand::MaxRelease {
            filename,
            output,
        } => releases::write_max(&filename, &output)?,
        Subcommand::Metadata {
            command,
        } => match command {
            MetadataCommand::Merge {
                basic,
                coordinates,
                rfam_hits,
                r2dt_hits,
                previous,
                output,
            } => {
                metadata::write_merge(
                    &basic,
                    &coordinates,
                    &rfam_hits,
                    &r2dt_hits,
                    &previous,
                    &output,
                )?;
            },
            MetadataCommand::Group {
                data_type,
                path,
                max_count,
                output,
            } => match data_type {
                Groupable::Basic => metadata::basic::group(&path, max_count, &output)?,
                Groupable::Coordinates => metadata::coordinate::group(&path, max_count, &output)?,
                Groupable::Previous => metadata::previous::group(&path, max_count, &output)?,
                Groupable::R2dtHits => metadata::r2dt_hit::group(&path, max_count, &output)?,
                Groupable::RfamHits => metadata::rfam_hit::group(&path, max_count, &output)?,
            },
        },
        Subcommand::GroupAccessions {
            path,
            min_index,
            max_index,
            output,
        } => accessions::group(&path, min_index, max_index, &output)?,
        Subcommand::Normalize {
            accessions,
            metadata,
            output,
        } => normalize::write(&accessions, &metadata, &output)?,
        Subcommand::Select {
            xrefs,
            known,
            output,
        } => releases::select_new(&xrefs, &known, &output)?,
    }

    Ok(())
}
