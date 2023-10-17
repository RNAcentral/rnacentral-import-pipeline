use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EditingEvent {
    pub id: usize,
    urs_taxid: String,
    has_editing_event: bool,
    chromosome: String,
    genomic_location: usize,
    repeat_type: String,
}

impl grouper::HasIndex for EditingEvent {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<EditingEvent>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}

impl EditingEvent {
    pub fn has_editing_event(&self) -> bool {
        self.has_editing_event
    }

    pub fn urs_taxid(&self) -> &str {
        &self.urs_taxid
    }

    pub fn genomic_location(&self) -> usize {
        self.genomic_location
    }

    pub fn repeat_type(&self) -> &str {
        &self.repeat_type
    }

    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }
}
