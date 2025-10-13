use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GoFlowLLMAnnotation {
    pub id: usize,
    urs_taxid: String,
    should_show_goflow: bool,
}

impl grouper::HasIndex for GoFlowLLMAnnotation {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<GoFlowLLMAnnotation>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}

impl GoFlowLLMAnnotation {
    pub fn should_show_goflow(&self) -> bool {
        self.should_show_goflow
    }
    pub fn urs_taxid(&self) -> &str {
        &self.urs_taxid
    }
}
