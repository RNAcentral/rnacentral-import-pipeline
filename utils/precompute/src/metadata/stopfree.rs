use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use rnc_core::grouper;

#[derive(Debug, Deserialize, Serialize)]
pub struct Stopfree {
    pub id: usize,
    pub is_protein_coding: Option<bool>,
}

impl grouper::HasIndex for Stopfree {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Stopfree>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}
