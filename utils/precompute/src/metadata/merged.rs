use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use crate::metadata::{
    basic::Basic,
    coordinate::Coordinate,
    previous::Previous,
    r2dt_hit::R2dtHit,
    rfam_hit::RfamHit,
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
}

impl Metadata {
    pub fn new(
        basic: Basic,
        coordinates: Vec<Coordinate>,
        rfam_hits: Vec<RfamHit>,
        r2dt_hits: Option<R2dtHit>,
        previous: Option<Previous>,
    ) -> Result<Self> {
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
        });
    }
}
