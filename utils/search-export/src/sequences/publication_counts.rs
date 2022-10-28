use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PublicationCount {
    pub id: usize,
    publication_count: usize,
}

impl grouper::HasIndex for PublicationCount {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<PublicationCount>(grouper::Criteria::ZeroOrOne, &path, 1, max, &output)
}

impl PublicationCount {
    pub fn publication_count(&self) -> usize {
        self.publication_count
    }
}
