use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use rnc_core::grouper;

#[derive(Debug, Deserialize, Serialize)]
pub struct Repeatmasker {
    pub id: usize,
    pub has_repeats: Option<bool>,
}

impl grouper::HasIndex for Repeatmasker {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Repeatmasker>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}
