use std::path::PathBuf;
use structopt::StructOpt;

use anyhow::Result;

use strum_macros::{
    Display,
    EnumIter,
    EnumString,
};

pub mod genes;
pub mod search_xml;
pub mod sequences;
pub mod utils;

#[derive(Debug, Display, PartialEq, Eq, Hash, EnumString, EnumIter)]
#[strum(ascii_case_insensitive, serialize_all = "kebab-case")]
pub enum Groupable {
    Base,
    Crs,
    Feedback,
    GoAnnotations,
    InteractingProteins,
    InteractingRnas,
    RegionInfo,
    Orfs,
    Precompute,
    QaStatus,
    R2dt,
    RfamHits,
    SoInfo,
}

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum GenesCommand {
    /// This command does a join between the sequence information, locus information and
    /// gene information. Any sequence not part of a gene is ignored. It then splits
    /// the resulting merged entries by assembly so it produces a directory of files,
    /// one per assmebly. These can then be merge together to produce complete gene
    /// entries. This assumes the locus and sequence information is sorted by a common
    /// id.
    SelectAndSplit {
        /// Filename to read the locus level information from
        #[structopt(parse(from_os_str))]
        locus: PathBuf,

        /// Filename to read the selected members from
        #[structopt(parse(from_os_str))]
        sequences: PathBuf,

        /// Output directory to write to, does not need to exist
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },

    /// This will take all gene member information and merge them into genes. This needs
    /// enough memory to fit information on all members from an assembly into memory
    /// at once.
    MergeAssembly {
        /// A file of all gene members for a given assembly. The file does not need to be
        /// sorted.
        #[structopt(parse(from_os_str))]
        members: PathBuf,

        /// Filename to write complete gene to
        #[structopt(parse(from_os_str))]
        output: PathBuf,
    },

    AsXml {
        /// File of the complete genes to read
        #[structopt(parse(from_os_str))]
        genes: PathBuf,

        /// Filename to write the xml data to
        #[structopt(parse(from_os_str))]
        xml_output: PathBuf,

        /// Filename to write the counts to
        #[structopt(parse(from_os_str))]
        count_output: PathBuf,
    },
}

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum SequenceCommand {
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
        orfs: PathBuf,

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

    /// Set of commands dealing with building data for sequence data export
    Sequences {
        #[structopt(subcommand)]
        command: SequenceCommand,
    },

    /// Commands dealing with building data for gene export
    Genes {
        #[structopt(subcommand)]
        command: GenesCommand,
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
            Groupable::Base => sequences::basic::group(&path, max_count, &output)?,
            Groupable::Crs => sequences::crs::group(&path, max_count, &output)?,
            Groupable::Feedback => sequences::feedback::group(&path, max_count, &output)?,
            Groupable::GoAnnotations => sequences::go_annotation::group(&path, max_count, &output)?,
            Groupable::InteractingProteins => {
                sequences::interacting_protein::group(&path, max_count, &output)?
            },
            Groupable::InteractingRnas => {
                sequences::interacting_rna::group(&path, max_count, &output)?
            },
            Groupable::RegionInfo => genes::region::group(&path, max_count, &output)?,
            Groupable::Precompute => sequences::precompute::group(&path, max_count, &output)?,
            Groupable::QaStatus => sequences::qa_status::group(&path, max_count, &output)?,
            Groupable::R2dt=> sequences::r2dt::group(&path, max_count, &output)?,
            Groupable::RfamHits => sequences::rfam_hit::group(&path, max_count, &output)?,
            Groupable::Orfs => sequences::orf::group(&path, max_count, &output)?,
            Groupable::SoInfo => Err(anyhow::anyhow!("May not group so info"))?,
        },
        Subcommand::Sequences {
            command,
        } => match command {
            SequenceCommand::Merge {
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
                orfs,
                so_term_tree,
                output,
            } => sequences::writers::write_merge(
                vec![
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
                    orfs,
                    so_term_tree,
                ],
                &output,
            )?,
            SequenceCommand::Normalize {
                accessions,
                metadata,
                output,
            } => sequences::writers::write(&accessions, &metadata, &output)?,
        },
        Subcommand::Genes {
            command,
        } => match command {
            GenesCommand::SelectAndSplit {
                locus,
                sequences,
                output,
            } => genes::writers::write_selected_members(&locus, &sequences, &output)?,
            GenesCommand::MergeAssembly {
                members,
                output,
            } => genes::writers::write_merged_members(&members, &output)?,
            GenesCommand::AsXml {
                genes,
                xml_output,
                count_output,
            } => genes::writers::write_search_files(&genes, &xml_output, &count_output)?,
        },
    }

    Ok(())
}
