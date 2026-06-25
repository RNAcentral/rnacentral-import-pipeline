use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

/// Presence of RepeatMasker-detected repetitive regions for a sequence. Only
/// urs_taxids with at least one repeat are emitted by the query, so this is a
/// zero-or-one (Optional) group used purely as a has_repetitive_region flag.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Repeatmasker {
    pub id: usize,
}

impl grouper::HasIndex for Repeatmasker {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Repeatmasker>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}
