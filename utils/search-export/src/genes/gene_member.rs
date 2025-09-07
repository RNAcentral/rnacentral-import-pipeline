use serde::{Deserialize, Serialize};

use crate::sequences::normalized::Normalized;

use super::region::UrsRegion;

#[derive(Debug, Deserialize, Serialize)]
pub struct GeneMember {
    region: UrsRegion,
    sequence: Normalized,
}

impl GeneMember {
    pub fn new(region: UrsRegion, sequence: Normalized) -> Self {
        assert!(
            region.id() == sequence.id(),
            "Region and sequence do not match ids, region: {} vs sequence {}",
            region.id(),
            sequence.id()
        );
        assert!(
            region.urs_taxid() == sequence.urs_taxid(),
            "Region and sequence do not match urs_taxids, region: {} vs sequence: {}",
            region.urs_taxid(),
            sequence.urs_taxid(),
        );

        Self {
            region,
            sequence,
        }
    }

    pub fn into_inner(self) -> (UrsRegion, Normalized) {
        (self.region, self.sequence)
    }

    pub fn gene_id(&self) -> usize {
        self.region.gene_id()
    }
}
