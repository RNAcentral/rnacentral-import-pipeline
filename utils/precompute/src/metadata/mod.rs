pub mod basic;
pub mod coordinate;
pub mod file_types;
pub mod merged;
pub mod orf;
pub mod previous;
pub mod r2dt_hit;
pub mod rfam_hit;
pub mod stopfree;
pub mod tcode;

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
    psql::PsqlJsonIterator,
};

pub use basic::Basic;
pub use coordinate::Coordinate;
pub use merged::Metadata;
pub use previous::Previous;
pub use r2dt_hit::R2dtHit;
pub use rfam_hit::RfamHit;
pub use stopfree::Stopfree;
pub use tcode::Tcode;

pub fn write_merge(
    basic_file: &Path,
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    orf_file: &Path,
    stopfree_file: &Path,
    tcode_file: &Path,
    output: &Path,
) -> Result<()> {
    let mut basic = PsqlJsonIterator::from_path(basic_file)?;
    let mut coordinates = PsqlJsonIterator::from_path(coordinate_file)?;
    let mut rfam_hits = PsqlJsonIterator::from_path(rfam_hits_file)?;
    let mut r2dt_hits = PsqlJsonIterator::from_path(r2dt_hits_file)?;
    let mut previous = PsqlJsonIterator::from_path(previous_file)?;
    let mut orfs = PsqlJsonIterator::from_path(orf_file)?;
    let mut stopfree = PsqlJsonIterator::from_path(stopfree_file)?;
    let mut tcode = PsqlJsonIterator::from_path(tcode_file)?;

    let mut output = rnc_utils::buf_writer(output)?;
    loop {
        match (
            basic.next(),
            coordinates.next(),
            rfam_hits.next(),
            r2dt_hits.next(),
            previous.next(),
            orfs.next(),
            stopfree.next(),
            tcode.next(),
        ) {
            (None, None, None, None, None, None, None, None) => break,
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
                Some(Optional {
                    id: id7,
                    data: stopfree,
                }),
                Some(Optional {
                    id: id8,
                    data: tcode,
                }),
            ) => {
                assert!(
                    id1 == id2
                        && id1 == id3
                        && id1 == id4
                        && id1 == id5
                        && id1 == id6
                        && id1 == id7
                        && id1 == id8,
                    "The data ids are out of sync at {}",
                    id1
                );

                let merged = Metadata::new(
                    basic,
                    coords,
                    rfam_hits,
                    r2dt_hit,
                    previous,
                    orfs,
                    stopfree,
                    tcode,
                )?;
                serde_json::to_writer(&mut output, &merged)?;
                writeln!(&mut output)?;
            },
            _ => return Err(anyhow!("Incorrect data format")),
        }
    }

    Ok(())
}
