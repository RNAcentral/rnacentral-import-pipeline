use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct R2dt {
    pub id: usize,
    secondary_structure_model: String,
    secondary_structure_source: String,
}

#[derive(Debug, PartialEq, Eq, Serialize)]
pub struct SecondaryStructure {
    model_name: String,
    model_source: String,
}

impl grouper::HasIndex for R2dt {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<R2dt>(grouper::Criteria::ZeroOrOne, &path, max, &output)
}
