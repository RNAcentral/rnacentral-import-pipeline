use std::path::PathBuf;
use std::str::FromStr;
use structopt::StructOpt;

use anyhow::Result;

pub mod normalize;

#[derive(Debug)]
enum Groupable {
    Basic,
    Crs,
    Feedback,
    GoAnnotations,
    InteractingProteins,
    InteractingRnas,
    QaStatus,
    R2dtHits,
    RfamHits,
}

impl FromStr for Groupable {
    type Err = String;

    fn from_str(raw: &str) -> Result<Self, Self::Err> {
        match raw {
            "basic" => Ok(Self::Basic),
            "crs" => Ok(Self::Crs),
            "feedback" => Ok(Self::Feedback),
            "go-annotations" => Ok(Self::GoAnnotations),
            "interacting-proteins" => Ok(Self::InteractingProteins),
            "interacting-rnas" => Ok(Self::InteractingRnas),
            "qa-status" => Ok(Self::QaStatus),
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
enum Subcommand {
    /// Merge each distinct type of metadata into single merged entries. This assumes that
    /// all files are sorted in the same way. Additionally, the xref file must have at
    /// least one entry for all possible urs_taxids.
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
            Subcommand::Group {
                data_type,
                path,
                max_count,
                output,
            } => match data_type {
                Groupable::Basic => normalize::basic::group(&path, max_count, &output)?,
                Groupable::Crs => normalize::crs::group(&path, max_count, &output)?,
                Groupable::Feedback => normalize::feedback::group(&path, max_count, &output)?,
                Groupable::GoAnnotations => normalize::go_annotation::group(&path, max_count, &output)?,
                Groupable::InteractingProteins => normalize::interacting_protein::group(&path, max_count, &output)?,
                Groupable::InteractingRnas => normalize::interacting_rna::group(&path, max_count, &output)?,
                Groupable::QaStatus => normalize::qa_status::group(&path, max_count, &output)?,
                Groupable::R2dtHits => normalize::r2dt::group(&path, max_count, &output)?,
                Groupable::RfamHits => normalize::rfam_hit::group(&path, max_count, &output)?,
            },
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
