use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use crate::metadata::{
    basic::Basic,
    coordinate::Coordinate,
    orf::{
        Orf,
        OrfInfo,
    },
    previous::Previous,
    r2dt_hit::R2dtHit,
    rfam_hit::RfamHit,
    stopfree::Stopfree,
    tcode::Tcode,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct Metadata {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    pub upi: String,
    pub taxid: usize,
    pub length: usize,
    pub coordinates: Vec<Coordinate>,
    pub previous: Option<Previous>,
    pub rfam_hits: Vec<RfamHit>,
    pub r2dt_hits: Option<R2dtHit>,
    pub orf_info: Option<OrfInfo>,
    pub possible_orf_stopfree: Option<bool>,
    pub possible_orf_tcode: Option<bool>,
}

impl Metadata {
    pub fn new(
        basic: Basic,
        coordinates: Vec<Coordinate>,
        rfam_hits: Vec<RfamHit>,
        r2dt_hits: Option<R2dtHit>,
        previous: Option<Previous>,
        orfs: Vec<Orf>,
        stopfree: Option<Stopfree>,
        tcode: Option<Tcode>,
    ) -> Result<Self> {
        if coordinates.len() > 0 {
            assert!(
                basic.urs_taxid == coordinates[0].urs_taxid,
                "Coordinates had incorrect urs_taxid {}, expected: {}",
                &coordinates[0].urs_taxid,
                &basic.urs_taxid
            )
        }

        if rfam_hits.len() > 0 {
            assert!(
                basic.urs_id == rfam_hits[0].urs_id,
                "Rfam Hits had incorrect urs_id {}, expected: {}",
                &rfam_hits[0].urs_id,
                &basic.urs_id
            )
        }

        let possible_orf_stopfree = stopfree.and_then(|s| s.is_protein_coding);
        let possible_orf_tcode = tcode.and_then(|t| t.is_protein_coding);

        return Ok(Self {
            id: basic.id,
            urs_id: basic.urs_id,
            urs_taxid: basic.urs_taxid.to_owned(),
            upi: basic.urs.to_owned(),
            taxid: basic.taxid,
            length: basic.length,
            coordinates,
            previous,
            rfam_hits,
            r2dt_hits,
            orf_info: orfs.into_iter().collect(),
            possible_orf_stopfree,
            possible_orf_tcode,
        });
    }
}
