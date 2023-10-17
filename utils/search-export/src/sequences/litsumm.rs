use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LitsummSummaries {
    pub id: usize,
    urs_taxid: String,
    should_show_litsumm: bool,
}

impl grouper::HasIndex for LitsummSummaries {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<LitsummSummaries>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}

impl LitsummSummaries {
    pub fn should_show_litsumm(&self) -> bool {
        self.should_show_litsumm
    }
    pub fn urs_taxid(&self) -> &str {
        &self.urs_taxid
    }
}
