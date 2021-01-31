pub mod chunks;
pub mod coordinate;
pub mod merged;
pub mod previous;
pub mod r2dt_hit;
pub mod rfam_hit;
pub mod xref;

use std::{
    fs::File,
    io::{
        prelude::*,
        BufWriter,
    },
    path::{
        Path,
        PathBuf,
    },
};

use anyhow::Result;

use itertools::Itertools;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

pub use coordinate::Coordinate;
pub use merged::Metadata;
pub use previous::Previous;
pub use r2dt_hit::R2dtHit;
pub use rfam_hit::RfamHit;
pub use xref::Xref;

use rnc_core::psql::JsonlIterator;

pub fn write_merge(
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    xref_file: &Path,
    output: &Path,
) -> Result<()> {
    let xrefs = JsonlIterator::from_path(xref_file)?;
    let xrefs = xrefs.group_by(|x: &Xref| x.urs_taxid.to_owned());
    let xrefs = xrefs.into_iter().assume_sorted_by_key();

    let coordinates = JsonlIterator::from_path(coordinate_file)?;
    let coordinates = coordinates.group_by(|c: &Coordinate| c.urs_taxid.to_owned());
    let coordinates = coordinates.into_iter().assume_sorted_by_key();

    let rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let rfam_hits = rfam_hits.group_by(|h: &RfamHit| h.urs_taxid.to_owned());
    let rfam_hits = rfam_hits.into_iter().assume_sorted_by_key();

    let r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let r2dt_hits = r2dt_hits.group_by(|h: &R2dtHit| h.urs_taxid.to_owned());
    let r2dt_hits = r2dt_hits.into_iter().assume_sorted_by_key();

    let previous = JsonlIterator::from_path(previous_file)?;
    let previous = previous.group_by(|p: &Previous| p.urs_taxid.to_owned());
    let previous = previous.into_iter().assume_sorted_by_key();

    let partial =
        xrefs.left_join(coordinates).left_join(rfam_hits).left_join(r2dt_hits).left_join(previous);

    let mut output = rnc_utils::buf_writer(output)?;
    for (id, data) in partial {
        let ((((xrefs, coordinates), rfam_hits), r2dt_hits), previous) = data;
        let norm = Metadata::new(id, xrefs, coordinates, rfam_hits, r2dt_hits, previous)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    Ok(())
}

pub fn write_splits(filename: &Path, chunks: &Path, output: &Path) -> Result<()> {
    let metadata = JsonlIterator::from_path(&filename)?;
    let chunks = chunks::load(chunks)?;
    let grouped = metadata.group_by(|metadata: &Metadata| {
        chunks
            .filename_of(metadata.ordering_index)
            .ok_or(format!("Could not find filename for {:?}", &metadata))
            .unwrap()
    });

    for (filename, items) in &grouped {
        let mut path = PathBuf::from(output);
        path.push(filename);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        for item in items {
            serde_json::to_writer(&mut writer, &item)?;
            writeln!(&mut writer)?;
        }
    }

    Ok(())
}
