use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;

use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Basic {
    pub id: usize,
    pub length: usize,
    pub md5: String,
    pub urs_taxid: String,
}

impl grouper::HasIndex for Basic {
    fn index(&self) -> usize {
        self.id
    }
}

impl Basic {
    pub fn urs_taxid(&self) -> &str {
        &self.urs_taxid
    }

    /// Get a mutable reference to the basic's id.
    pub fn id_mut(&mut self) -> &mut usize {
        &mut self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Basic>(grouper::Criteria::ExactlyOne, &path, 1, max, &output)
}
