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
pub struct InteractingRna {
    pub id: usize,
    interacting_rna_id: String,
    urs: Option<String>,
    methods: Vec<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct InteractingRnaVec {
    interacting_rna_id: HashSet<String>,
    urs: HashSet<String>,
    methods: HashSet<String>,
}

impl Default for InteractingRnaVec {
    fn default() -> Self {
        Self {
            interacting_rna_id: HashSet::new(),
            urs: HashSet::new(),
            methods: HashSet::new(),
        }
    }
}

impl FromIterator<InteractingRna> for InteractingRnaVec {
    fn from_iter<I: IntoIterator<Item = InteractingRna>>(iter: I) -> Self {
        let mut value = InteractingRnaVec::default();

        for i in iter {
            value.interacting_rna_id.insert(i.interacting_rna_id);
            if i.urs.is_some() {
                value.urs.insert(i.urs.unwrap());
            }
            value.methods.extend(i.methods);
        }

        value
    }
}

impl grouper::HasIndex for InteractingRna {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<InteractingRna>(grouper::Criteria::AnyNumber, &path, max, &output)
}
