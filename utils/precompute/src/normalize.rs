use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use itertools::Itertools;

use anyhow::Result;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use crate::{
    accessions::{
        Accession,
        RawAccessionEntry,
    },
    metadata::{
        coordinate::Coordinate,
        merged::Metadata,
        previous::Previous,
        r2dt_hit::R2dtHit,
        rfam_hit::RfamHit,
    },
};
use rnc_core::psql::JsonlIterator;

#[derive(Debug, Deserialize, Serialize)]
pub struct Normalized {
    upi: String,
    taxid: usize,
    length: usize,
    last_release: usize,
    coordinates: Vec<Coordinate>,
    accessions: Vec<Accession>,
    deleted: bool,
    previous: Option<Previous>,
    rfam_hits: Vec<RfamHit>,
    r2dt_hits: Vec<R2dtHit>,
}

impl Normalized {
    fn new(
        raw_accessions: impl Iterator<Item = RawAccessionEntry>,
        metadata: Metadata,
    ) -> Result<Self> {
        let accessions = raw_accessions
            .map(|a: RawAccessionEntry| {
                let is_active = metadata
                    .active_mapping
                    .get(&a.accession)
                    .ok_or(format!("Missing {:} in {:?}", &a.accession, &metadata.active_mapping))
                    .unwrap();
                Accession::from((a, is_active.clone()))
            })
            .collect();
        return Ok(Self {
            upi: metadata.upi,
            taxid: metadata.taxid,
            length: metadata.length,
            last_release: metadata.last_release,
            coordinates: metadata.coordinates,
            accessions,
            deleted: metadata.deleted,
            previous: metadata.previous,
            rfam_hits: metadata.rfam_hits,
            r2dt_hits: metadata.r2dt_hits,
        });
    }
}

pub fn write(accession_file: &Path, metadata_file: &Path, output: &Path) -> Result<()> {
    let accessions = JsonlIterator::from_path(accession_file)?;
    let accessions = accessions.group_by(|a: &RawAccessionEntry| (a.urs_id, a.id));
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let metadata = JsonlIterator::from_path(metadata_file)?;
    let metadata = metadata.map(|m: Metadata| ((m.urs_id, m.id), m));
    let metadata = metadata.into_iter().assume_sorted_by_key();

    let partial = accessions.join(metadata);

    let mut output = rnc_utils::buf_writer(output)?;
    for (_id, data) in partial {
        let (accessions, metadata) = data;
        let norm = Normalized::new(accessions, metadata)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    Ok(())
}
