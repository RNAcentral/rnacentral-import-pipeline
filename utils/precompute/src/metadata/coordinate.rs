use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use rnc_core::grouper;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Coordinate {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    assembly_id: String,
    chromosome: String,
    strand: i8,
    start: usize,
    stop: usize,
}

impl grouper::HasIndex for Coordinate {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Coordinate>(grouper::Criteria::AnyNumber, &path, max, &output)
}
