pub mod basic;
pub mod coordinate;
pub mod merged;
pub mod orf;
pub mod previous;
pub mod r2dt_hit;
pub mod rfam_hit;

use std::{
    io::Write,
    path::Path,
};

use anyhow::{
    anyhow,
    Result,
};

use rnc_core::{
    grouper::Grouped::{
        Multiple,
        Optional,
        Required,
    },
    psql::JsonlIterator,
};

pub use basic::Basic;
pub use coordinate::Coordinate;
pub use merged::Metadata;
pub use previous::Previous;
pub use r2dt_hit::R2dtHit;
pub use rfam_hit::RfamHit;

pub fn write_merge(
    basic_file: &Path,
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    orf_file: &Path,
    output: &Path,
) -> Result<()> {
    let mut basic = JsonlIterator::from_path(basic_file)?;
    let mut coordinates = JsonlIterator::from_path(coordinate_file)?;
    let mut rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let mut r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let mut previous = JsonlIterator::from_path(previous_file)?;
    let mut orfs = JsonlIterator::from_path(orf_file)?;

    let mut output = rnc_utils::buf_writer(output)?;
    loop {
        match (
            basic.next(),
            coordinates.next(),
            rfam_hits.next(),
            r2dt_hits.next(),
            previous.next(),
            orfs.next(),
        ) {
            (None, None, None, None, None, None) => break,
            (
                Some(Required {
                    id: id1,
                    data: basic,
                }),
                Some(Multiple {
                    id: id2,
                    data: coords,
                }),
                Some(Multiple {
                    id: id3,
                    data: rfam_hits,
                }),
                Some(Optional {
                    id: id4,
                    data: r2dt_hit,
                }),
                Some(Optional {
                    id: id5,
                    data: previous,
                }),
                Some(Multiple {
                    id: id6,
                    data: orfs,
                }),
            ) => {
                assert!(
                    id1 == id2 && id1 == id3 && id1 == id4 && id1 == id5 && id1 == id6,
                    "The data ids are out of sync at {}",
                    id1
                );

                let merged = Metadata::new(basic, coords, rfam_hits, r2dt_hit, previous, orfs)?;
                serde_json::to_writer(&mut output, &merged)?;
                writeln!(&mut output)?;
            },
            _ => return Err(anyhow!("Incorrect data format")),
        }
    }

    Ok(())
}
