use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct QaStatus {
    pub id: usize,
    has_issue: bool,
    possible_contamination: bool,
    incomplete_sequence: bool,
    missing_rfam_match: bool,
}

impl grouper::HasIndex for QaStatus {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<QaStatus>(grouper::Criteria::ExactlyOne, &path, 1, max, &output)
}
