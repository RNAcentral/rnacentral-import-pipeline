use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use anyhow::Result;
use rnc_core::grouper;

/// A disease a source database links to this sequence. The term and its name are
/// missing when the source's free text could not be mapped to an ontology term.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DiseaseAnnotation {
    pub id: usize,
    urs_taxid: String,
    disease_term: Option<String>,
    disease_name: Option<String>,
    disease_text: String,
    assigned_by: String,
}

impl grouper::HasIndex for DiseaseAnnotation {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<DiseaseAnnotation>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}
