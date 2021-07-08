use std::path::{
    Path,
    PathBuf,
};

use anyhow::{
    anyhow,
    Result,
};

use itertools::Itertools;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use rnc_core::psql::JsonlIterator;

use self::file_joiner::{
    FileJoinerBuilder,
    FileTypes,
};

pub mod accession;
pub mod basic;
pub mod crs;
pub mod entry;
pub mod feedback;
pub mod file_joiner;
pub mod go_annotation;
pub mod interacting_protein;
pub mod interacting_rna;
pub mod orf;
pub mod precompute;
pub mod qa_status;
pub mod r2dt;
pub mod rfam_hit;
pub mod so_tree;

pub fn write_merge(
    base_file: PathBuf,
    crs_file: PathBuf,
    feedback_file: PathBuf,
    go_annotations_file: PathBuf,
    interacting_proteins_file: PathBuf,
    interacting_rnas_file: PathBuf,
    precompute_file: PathBuf,
    qa_status_file: PathBuf,
    r2dt_hits_file: PathBuf,
    rfam_hits_file: PathBuf,
    orf_file: PathBuf,
    so_term_tree_file: PathBuf,
    output_file: &Path,
) -> Result<()> {
    let mut builder = FileJoinerBuilder::default();
    builder
        .file(FileTypes::Basic, base_file)
        .file(FileTypes::Crs, crs_file)
        .file(FileTypes::Feedback, feedback_file)
        .file(FileTypes::GoAnnotations, go_annotations_file)
        .file(FileTypes::InteractingProteins, interacting_proteins_file)
        .file(FileTypes::InteractingRnas, interacting_rnas_file)
        .file(FileTypes::Orfs, orf_file)
        .file(FileTypes::Precompute, precompute_file)
        .file(FileTypes::QaStatus, qa_status_file)
        .file(FileTypes::R2dtHits, r2dt_hits_file)
        .file(FileTypes::RfamHits, rfam_hits_file)
        .file(FileTypes::SoInfo, so_term_tree_file);

    let joiner = builder.build()?;
    let mut writer = rnc_utils::buf_writer(output_file)?;
    for raw in joiner {
        let raw = raw?;
        serde_json::to_writer(&mut writer, &raw)?;
        writeln!(&mut writer)?;
    }

    Ok(())
}

pub fn write(accession_file: &Path, metadata_file: &Path, output_file: &Path) -> Result<()> {
    let mut writer = rnc_utils::buf_writer(output_file)?;

    let metadata = JsonlIterator::from_path(metadata_file)?;
    let metadata = metadata.map(|r: entry::Raw| (r.id(), r));
    let metadata = metadata.into_iter().assume_sorted_by_key();

    let accessions = JsonlIterator::from_path(accession_file)?;
    let accessions = accessions.group_by(|a: &accession::RawAccession| a.id);
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let merged = accessions.join(metadata);

    let mut seen = false;
    for (_id, entry) in merged {
        let (accessions, metadata) = entry;
        let normalized = entry::Normalized::new(metadata, accessions.collect())?;
        serde_json::to_writer(&mut writer, &normalized)?;
        writeln!(&mut writer)?;
        seen = true;
    }

    if !seen {
        return Err(anyhow!("Failed to find any normalized data"));
    }

    Ok(())
}
