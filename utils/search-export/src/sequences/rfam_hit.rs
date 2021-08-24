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

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RfamHit {
    pub id: usize,
    rfam_ids: String,
    rfam_family_names: String,
    rfam_clans: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RfamHitVec {
    rfam_ids: HashSet<String>,
    rfam_family_names: HashSet<String>,
    rfam_clans: HashSet<String>,
}

impl Default for RfamHitVec {
    fn default() -> Self {
        Self {
            rfam_ids: HashSet::new(),
            rfam_family_names: HashSet::new(),
            rfam_clans: HashSet::new(),
        }
    }
}

impl FromIterator<RfamHit> for RfamHitVec {
    fn from_iter<I: IntoIterator<Item = RfamHit>>(iter: I) -> Self {
        let mut value = RfamHitVec::default();

        for i in iter {
            value.rfam_ids.insert(i.rfam_ids);
            value.rfam_family_names.insert(i.rfam_family_names);
            value.rfam_clans.extend(i.rfam_clans);
        }

        value
    }
}

impl grouper::HasIndex for RfamHit {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<RfamHit>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}
