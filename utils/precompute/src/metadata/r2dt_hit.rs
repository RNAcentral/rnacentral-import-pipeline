use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use rnc_core::grouper;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct R2dtHit {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    model_id: usize,
    model_name: String,
    model_source: String,
    model_so_term: Option<String>,
    sequence_coverage: Option<f64>,
    model_coverage: Option<f64>,
    sequence_basepairs: Option<usize>,
    model_basepairs: Option<usize>,
}

impl grouper::HasIndex for R2dtHit {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<R2dtHit>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}
