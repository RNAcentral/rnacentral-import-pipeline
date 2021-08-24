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
pub struct InteractingProtein {
    pub id: usize,
    interacting_protein_id: String,
    synonyms: Option<Vec<String>>,
    label: Option<String>,
    relationship: String,
    methods: Vec<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct InteractingProteinVec {
    interacting_protein_id: HashSet<String>,
    synonyms: HashSet<String>,
    label: HashSet<String>,
    relationship: HashSet<String>,
    methods: HashSet<String>,
}

impl Default for InteractingProteinVec {
    fn default() -> Self {
        Self {
            interacting_protein_id: HashSet::new(),
            synonyms: HashSet::new(),
            label: HashSet::new(),
            relationship: HashSet::new(),
            methods: HashSet::new(),
        }
    }
}

impl FromIterator<InteractingProtein> for InteractingProteinVec {
    fn from_iter<I: IntoIterator<Item = InteractingProtein>>(iter: I) -> Self {
        let mut value = InteractingProteinVec::default();

        for i in iter {
            value.interacting_protein_id.insert(i.interacting_protein_id);
            if i.synonyms.is_some() {
                value.synonyms.extend(i.synonyms.unwrap());
            }
            if i.label.is_some() {
                value.label.insert(i.label.unwrap());
            }
            value.relationship.insert(i.relationship);
            value.methods.extend(i.methods);
        }

        value
    }
}

impl grouper::HasIndex for InteractingProtein {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<InteractingProtein>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}
