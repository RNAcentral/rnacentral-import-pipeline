use std::{
    collections::HashSet,
    iter::FromIterator,
};

use serde::{
    Deserialize,
    Serialize,
};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Crs {
    pub id: usize,
    crs_id: String,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CrsVec {
    crs_id: HashSet<String>,
}

impl Default for CrsVec {
    fn default() -> Self {
        Self {
            crs_id: HashSet::new(),
        }
    }
}

impl FromIterator<Crs> for CrsVec {
    fn from_iter<I: IntoIterator<Item = Crs>>(iter: I) -> Self {
        let mut c = CrsVec::default();

        for i in iter {
            c.crs_id.insert(i.crs_id);
        }

        c
    }
}
