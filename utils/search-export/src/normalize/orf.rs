use std::{
    collections::HashSet,
    iter::FromIterator,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;
use rnc_core::grouper;

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Orf {
    pub id: usize,
    source: String,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrfVec {
    orf_sources: HashSet<String>,
}

impl Default for OrfVec {
    fn default() -> Self {
        Self {
            orf_sources: HashSet::default(),
        }
    }
}

impl FromIterator<Orf> for OrfVec {
    fn from_iter<I: IntoIterator<Item = Orf>>(orfs: I) -> Self {
        let mut summary = OrfVec::default();
        for orf in orfs {
            summary.orf_sources.insert(orf.source);
        }
        summary
    }
}

impl grouper::HasIndex for Orf {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Orf>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}
