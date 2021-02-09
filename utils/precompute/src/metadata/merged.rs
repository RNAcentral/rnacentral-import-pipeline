use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use crate::metadata::{
    coordinate::Coordinate,
    previous::Previous,
    r2dt_hit::R2dtHit,
    rfam_hit::RfamHit,
    basic::Basic,
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
    pub r2dt_hits: Vec<R2dtHit>,
}

impl Metadata {
    pub fn new(
        basic: Basic,
        raw_coordinates: Option<impl Iterator<Item = Coordinate>>,
        raw_rfam_hits: Option<impl Iterator<Item = RfamHit>>,
        raw_r2dt_hits: Option<impl Iterator<Item = R2dtHit>>,
        raw_previous: Option<impl Iterator<Item = Previous>>,
    ) -> Result<Self> {
        let coordinates = raw_coordinates.into_iter().flatten().collect();
        let rfam_hits = raw_rfam_hits.into_iter().flatten().collect();
        let r2dt_hits = raw_r2dt_hits.into_iter().flatten().collect();
        let previous: Vec<Previous> = raw_previous.into_iter().flatten().collect();
        let previous = match previous.len() {
            0 => None,
            _ => Some(previous[0].clone()),
        };

        return Ok(Self {
            id: basic.id,
            urs_id: basic.urs_id,
            urs_taxid: basic.urs_taxid,
            upi: basic.urs,
            taxid: basic.taxid,
            length: basic.length,
            coordinates,
            previous,
            rfam_hits,
            r2dt_hits,
        });
    }
}
