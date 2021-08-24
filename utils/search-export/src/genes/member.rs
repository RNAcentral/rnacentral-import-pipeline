use serde::{
    Deserialize,
    Serialize,
};

use crate::sequences::entry::Normalized;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct GeneMembers {
    locus_id: usize,
    sequences: Vec<Normalized>,
}

impl GeneMembers {
    pub fn new(sequences: Vec<Normalized>) -> Self {
        Self {
            locus_id: sequences[0].locus_id().unwrap(),
            sequences,
        }
    }

    pub fn locus_id(&self) -> usize {
        return self.locus_id;
    }

    pub fn sequences(&self) -> &Vec<Normalized> {
        &self.sequences
    }
}
