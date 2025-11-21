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

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Orf {
    pub id: usize,
    pub urs_id: usize,
    urs_taxid: String,
    source: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct OrfInfo {
    sources: HashSet<String>,
}

impl Default for OrfInfo {
    fn default() -> Self {
        Self {
            sources: HashSet::new(),
        }
    }
}

impl FromIterator<Orf> for Option<OrfInfo> {
    fn from_iter<I: IntoIterator<Item = Orf>>(orfs: I) -> Self {
        let mut info = OrfInfo::default();

        for orf in orfs {
            info.sources.insert(orf.source);
        }

        match info.sources.is_empty() {
            true => None,
            false => Some(info),
        }
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
