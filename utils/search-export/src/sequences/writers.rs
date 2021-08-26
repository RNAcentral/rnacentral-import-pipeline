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

use crate::sequences::{
    accession,
    normalized::Normalized,
    file_joiner::FileJoinerBuilder,
};

use super::raw::Raw;

pub fn write_merge(files: Vec<PathBuf>, output_file: &Path) -> Result<()> {
    let builder = FileJoinerBuilder::from_files(files)?;
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
    let metadata = metadata.map(|r: Raw| (r.id(), r));
    let metadata = metadata.into_iter().assume_sorted_by_key();

    let accessions = JsonlIterator::from_path(accession_file)?;
    let accessions = accessions.group_by(|a: &accession::RawAccession| a.id);
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let merged = accessions.join(metadata);

    let mut seen = false;
    for (_id, entry) in merged {
        let (accessions, metadata) = entry;
        let normalized = Normalized::new(metadata, accessions.collect())?;
        serde_json::to_writer(&mut writer, &normalized)?;
        writeln!(&mut writer)?;
        seen = true;
    }

    if !seen {
        return Err(anyhow!("Failed to find any normalized data"));
    }

    Ok(())
}
