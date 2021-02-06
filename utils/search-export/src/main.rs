use std::path::PathBuf;
use structopt::StructOpt;

use anyhow::Result;

pub mod normalize;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Subcommand {
    /// Merge each distinct type of metadata into single merged entries. This assumes that
    /// all files are sorted in the same way. Additionally, the xref file must have at
    /// least one entry for all possible urs_taxids.
    Merge {
        #[structopt(parse(from_os_str))]
        base: PathBuf,

        #[structopt(parse(from_os_str))]
        crs: PathBuf,

        #[structopt(parse(from_os_str))]
        feedback: PathBuf,

        #[structopt(parse(from_os_str))]
        go_annotations: PathBuf,

        #[structopt(parse(from_os_str))]
        interacting_proteins: PathBuf,

        #[structopt(parse(from_os_str))]
        interacting_rnas: PathBuf,

        #[structopt(parse(from_os_str))]
        precompute: PathBuf,

        #[structopt(parse(from_os_str))]
        qa_status: PathBuf,

        #[structopt(parse(from_os_str))]
        r2dt_hits: PathBuf,

        #[structopt(parse(from_os_str))]
        rfam_hits: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename of the SO term tree metadata.
        so_term_tree: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename to write the results to, '-' means stdout
        output: PathBuf,
    },

    Normalize {
        #[structopt(parse(from_os_str))]
        /// Filename of the accessions
        accessions: PathBuf,

        #[structopt(parse(from_os_str))]
        /// Filename for the merged metadata
        metadata: PathBuf,

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

fn main() -> Result<()> {
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
        Subcommand::Merge {
            base,
            crs,
            feedback,
            go_annotations,
            interacting_proteins,
            interacting_rnas,
            precompute,
            qa_status,
            r2dt_hits,
            rfam_hits,
            so_term_tree,
            output,
        } => normalize::write_merge(
            &base,
            &crs,
            &feedback,
            &go_annotations,
            &interacting_proteins,
            &interacting_rnas,
            &precompute,
            &qa_status,
            &r2dt_hits,
            &rfam_hits,
            &so_term_tree,
            &output,
        )?,
        Subcommand::Normalize {
            accessions,
            metadata,
            output,
        } => normalize::write(&accessions, &metadata, &output)?,
    }

    Ok(())
}
