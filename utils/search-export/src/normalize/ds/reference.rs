use std::{
    collections::HashSet,
    iter::FromIterator,
};

use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Reference {
    pub id: usize,
    authors: Option<String>,
    journal: Option<String>,
    pub_title: Option<String>,
    pub_id: u64,
    pubmed_id: Option<String>,
    doi: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReferenceVec {
    authors: HashSet<String>,
    journals: HashSet<String>,
    pub_titles: HashSet<String>,
    pub_ids: HashSet<u64>,
    pubmed_ids: HashSet<String>,
    dois: HashSet<String>,
}

impl Default for ReferenceVec {
    fn default() -> Self {
        Self {
            authors: HashSet::new(),
            journals: HashSet::new(),
            pub_titles: HashSet::new(),
            pub_ids: HashSet::new(),
            pubmed_ids: HashSet::new(),
            dois: HashSet::new(),
        }
    }
}

impl FromIterator<Reference> for ReferenceVec {
    fn from_iter<I: IntoIterator<Item = Reference>>(iter: I) -> Self {
        let mut value = ReferenceVec::default();

        for i in iter {
            if i.authors.is_some() {
                let authors = i.authors.unwrap();
                let authors = authors.split(", ").filter(|a| !a.is_empty()).map(|s| s.to_string());
                value.authors.extend(authors);
            }
            if i.journal.is_some() {
                value.journals.insert(i.journal.unwrap());
            }
            if i.pub_title.is_some() {
                value.pub_titles.insert(i.pub_title.unwrap());
            }
            value.pub_ids.insert(i.pub_id);
            if i.pubmed_id.is_some() {
                value.pubmed_ids.insert(i.pubmed_id.unwrap());
            }
            if i.doi.is_some() {
                value.dois.insert(i.doi.unwrap());
            }
        }

        value
    }
}
