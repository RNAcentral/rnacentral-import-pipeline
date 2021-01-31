use std::path::PathBuf;
use structopt::StructOpt;

pub mod accessions;
pub mod metadata;
pub mod normalize;
pub mod releases;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum MetadataCommand {
    /// Merge each distinct type of metadata into single merged entries. THis assumes that
    /// all files are sorted in the same way. Additionally, the xref file must have at
    /// least one entry for all possible urs_taxids.
    Merge {
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
        /// Filename of the raw xref file
        xref: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    /// Split the metadata into a series of chunks and write the chunks into separate
    /// files. Doing this makes the normalzing step faster if you use the same chunks
    /// as in accession chunks.
    Chunk {
        #[structopt(parse(from_os_str))]
        /// The name of the metadata to read.
        metadata: PathBuf,

        #[structopt(parse(from_os_str))]
        /// The name of the metadata to read.
        chunk_file: PathBuf,

        #[structopt(parse(from_os_str), default_value = ".")]
        /// Base directory to write to.
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
                coordinates,
                rfam_hits,
                r2dt_hits,
                previous,
                xref,
                output,
            } => {
                metadata::write_merge(
                    &coordinates,
                    &rfam_hits,
                    &r2dt_hits,
                    &previous,
                    &xref,
                    &output,
                )?;
            },
            MetadataCommand::Chunk {
                metadata,
                chunk_file,
                output,
            } => metadata::write_splits(&metadata, &chunk_file, &output)?,
        },
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
