pub mod coordinate;
pub mod merged;
pub mod previous;
pub mod r2dt_hit;
pub mod rfam_hit;
pub mod basic;

use std::{
    io::Write,
    path::Path,
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
pub use basic::Basic;

use rnc_core::psql::JsonlIterator;

pub fn write_merge(
    basic_file: &Path,
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    output: &Path,
) -> Result<()> {
    let basic = JsonlIterator::from_path(basic_file)?;
    let basic = basic.map(|b: Basic| ((b.urs_id, b.id), b));
    let basic = basic.into_iter().assume_sorted_by_key();

    let coordinates = JsonlIterator::from_path(coordinate_file)?;
    let coordinates = coordinates.group_by(|c: &Coordinate| (c.urs_id, c.id));
    let coordinates = coordinates.into_iter().assume_sorted_by_key();

    let rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let rfam_hits = rfam_hits.group_by(|h: &RfamHit| (h.urs_id, h.id));
    let rfam_hits = rfam_hits.into_iter().assume_sorted_by_key();

    let r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let r2dt_hits = r2dt_hits.group_by(|h: &R2dtHit| (h.urs_id, h.id));
    let r2dt_hits = r2dt_hits.into_iter().assume_sorted_by_key();

    let previous = JsonlIterator::from_path(previous_file)?;
    let previous = previous.group_by(|p: &Previous| (p.urs_id, p.id));
    let previous = previous.into_iter().assume_sorted_by_key();

    let partial =
        basic.left_join(coordinates).left_join(rfam_hits).left_join(r2dt_hits).left_join(previous);

    let mut output = rnc_utils::buf_writer(output)?;
    for (_ids, data) in partial {
        let ((((basic, coordinates), rfam_hits), r2dt_hits), previous) = data;
        let norm = Metadata::new(basic, coordinates, rfam_hits, r2dt_hits, previous)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    Ok(())
}
