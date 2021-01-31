use serde::{
    Deserialize,
    Serialize,
};

use std::collections::HashMap;

use anyhow::Result;

use crate::metadata::{
    coordinate::Coordinate,
    previous::Previous,
    r2dt_hit::R2dtHit,
    rfam_hit::RfamHit,
    xref::Xref,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct Metadata {
    pub urs_taxid: String,
    pub upi: String,
    pub taxid: usize,
    pub length: usize,
    pub last_release: usize,
    pub coordinates: Vec<Coordinate>,
    pub deleted: bool,
    pub previous: Option<Previous>,
    pub rfam_hits: Vec<RfamHit>,
    pub r2dt_hits: Vec<R2dtHit>,
    pub ordering_index: usize,
    pub active_mapping: HashMap<String, bool>,
}

impl Metadata {
    pub fn new(
        urs_taxid: String,
        raw_xrefs: impl Iterator<Item = Xref>,
        raw_coordinates: Option<impl Iterator<Item = Coordinate>>,
        raw_rfam_hits: Option<impl Iterator<Item = RfamHit>>,
        raw_r2dt_hits: Option<impl Iterator<Item = R2dtHit>>,
        raw_previous: Option<impl Iterator<Item = Previous>>,
    ) -> Result<Self> {
        let parts: Vec<&str> = urs_taxid.split("_").collect();
        let upi = parts[0].to_owned();
        let taxid = parts[1].parse::<usize>().unwrap();

        let xrefs: Vec<Xref> = raw_xrefs.collect();
        let ordering_index = xrefs.get(0).unwrap().ordering_index;
        let mut active_mapping: HashMap<String, bool> = HashMap::new();
        for xref in &xrefs {
            active_mapping.insert(xref.accession.to_owned(), xref.is_active());
        }

        let coordinates = raw_coordinates.into_iter().flatten().collect();
        let rfam_hits = raw_rfam_hits.into_iter().flatten().collect();
        let r2dt_hits = raw_r2dt_hits.into_iter().flatten().collect();
        let previous: Vec<Previous> = raw_previous.into_iter().flatten().collect();
        let previous = match previous.len() {
            0 => None,
            _ => Some(previous[0].clone()),
        };

        let length = xrefs[0].length;
        let last_release = xrefs.iter().map(|x| x.last_release).max().unwrap();

        let deleted = xrefs.iter().all(|x| x.deleted == "Y");

        return Ok(Self {
            urs_taxid,
            upi,
            taxid,
            length,
            last_release,
            coordinates,
            deleted,
            previous,
            rfam_hits,
            r2dt_hits,
            ordering_index,
            active_mapping,
        });
    }
}
