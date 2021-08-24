use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LocusInfo {
    pub id: usize,
    locus_id: usize,
    membership_status: String,
}

impl grouper::HasIndex for LocusInfo {
    fn index(&self) -> usize {
        self.id
    }
}

impl LocusInfo {
    pub fn locus_id(&self) -> usize {
        self.locus_id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<LocusInfo>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}
