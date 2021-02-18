use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use crate::metadata::grouper;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct RfamHit {
    pub id: usize,
    pub urs_id: usize,
    urs_taxid: String,
    rfam_hit_id: usize,
    model: String,
    model_rna_type: Option<String>,
    model_domain: Option<String>,
    model_name: String,
    model_long_name: String,
    model_completeness: f64,
    model_start: usize,
    model_stop: usize,
    sequence_completeness: f64,
    sequence_start: usize,
    sequence_stop: usize,
}

impl grouper::HasIndex for RfamHit {
    fn index(&self) -> usize {
        self.urs_id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<RfamHit>(grouper::Criteria::AnyNumber, &path, max, &output)
}
