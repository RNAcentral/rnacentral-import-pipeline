use serde::{
    Deserialize,
    Serialize,
};

use crate::sequences::{
    accession::CrossReference,
    normalized::Normalized,
};

use super::region::Region;

#[derive(Debug, Deserialize, Serialize)]
pub struct GeneMember {
    region: Region,
    sequence: Normalized,
}

impl GeneMember {
    pub fn new(region: Region, sequence: Normalized) -> Self {
        assert!(region.id() == sequence.id(), "Region and sequence do not match ids");
        assert!(
            region.urs_taxid() == sequence.urs_taxid(),
            "Region and sequence do not match sequences"
        );

        Self {
            region,
            sequence,
        }
    }

    pub fn locus_id(&self) -> usize {
        self.region.locus_id()
    }

    pub fn name(&self) -> &str {
        self.region.name()
    }

    pub fn member_count(&self) -> &usize {
        &self.region.member_count()
    }

    pub fn assembly_id(&self) -> &str {
        self.region.assembly_id()
    }

    pub fn cross_references(&self) -> &[CrossReference] {
        self.sequence.cross_references()
    }

    pub fn urs_taxid(&self) -> &str {
        self.sequence.urs_taxid()
    }

    pub fn description(&self) -> &str {
        self.sequence.description()
    }
}
