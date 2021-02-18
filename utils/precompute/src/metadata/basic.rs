use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use crate::metadata::grouper;

#[derive(Debug, Deserialize, Serialize)]
pub struct Basic {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    pub urs: String,
    pub taxid: usize,
    pub length: usize,
}

impl grouper::HasIndex for Basic {
    fn index(&self) -> usize {
        self.urs_id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Basic>(grouper::Criteria::ExactlyOne, &path, max, &output)
}
