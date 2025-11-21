use serde::{
    Deserialize,
    Serialize,
};

use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoAnnotation {
    pub id: usize,
    go_term_id: String,
    qualifier: String,
    go_name: String,
    assigned_by: String,
}

impl grouper::HasIndex for GoAnnotation {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<GoAnnotation>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}
