use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use rnc_core::grouper;

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct Previous {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    upi: String,
    taxid: usize,
    databases: Option<String>,
    description: Option<String>,
    has_coordinates: Option<bool>,
    is_active: Option<bool>,
    last_release: Option<usize>,
    rna_type: Option<String>,
    short_description: Option<String>,
    so_rna_type: Option<String>,
}

impl grouper::HasIndex for Previous {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Previous>(grouper::Criteria::ZeroOrOne, &path, max, &output)
}
