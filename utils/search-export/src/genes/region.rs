use std::{
    collections::HashMap,
    fs::File,
    iter::FromIterator,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;
use rnc_core::{
    grouper::{
        self,
        Grouped,
    },
    psql::JsonlIterator,
};

use crate::utils::set_or_check;

#[derive(Debug, Deserialize, Serialize)]
pub struct SequenceWithRegions {
    id: usize,
    urs_taxid: String,
    regions: HashMap<String, Vec<Region>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Region {
    id: usize,
    urs_taxid: String,
    locus_id: usize,
    name: String,
    assembly_id: String,
    member_count: usize,
    membership_status: String,
}

impl Region {
    /// Get a reference to the region's locus id.
    pub fn locus_id(&self) -> usize {
        self.locus_id
    }

    /// Get a reference to the region's assembly id.
    pub fn assembly_id(&self) -> &str {
        self.assembly_id.as_str()
    }

    /// Get a reference to the region's name.
    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    /// Get a reference to the region's member count.
    pub fn member_count(&self) -> &usize {
        &self.member_count
    }

    /// Get a reference to the region's id.
    pub fn id(&self) -> &usize {
        &self.id
    }

    /// Get a reference to the region's urs taxid.
    pub fn urs_taxid(&self) -> &str {
        self.urs_taxid.as_str()
    }
}

pub struct RegionGrouper {
    iter: JsonlIterator<File, Grouped<Region>>,
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Region>(grouper::Criteria::AnyNumber, &path, 1, max, &output)
}

impl grouper::HasIndex for Region {
    fn index(&self) -> usize {
        self.id
    }
}

impl SequenceWithRegions {
    /// Get a reference to the sequence locus's id.
    pub fn id(&self) -> &usize {
        &self.id
    }

    /// Get a reference to the sequence with regions's regions.
    pub fn regions_by_assembly(&self) -> &HashMap<String, Vec<Region>> {
        &self.regions
    }
}

impl RegionGrouper {
    pub fn new(iter: JsonlIterator<File, Grouped<Region>>) -> Self {
        Self {
            iter,
        }
    }
}

impl FromIterator<Region> for SequenceWithRegions {
    fn from_iter<T: IntoIterator<Item = Region>>(iter: T) -> Self {
        let mut regions: HashMap<String, Vec<Region>> = HashMap::new();
        let mut id = None;
        let mut urs_taxid = None;
        for item in iter {
            set_or_check(&mut id, item.id);
            set_or_check(&mut urs_taxid, item.urs_taxid.to_string());
            let current = regions.entry(item.assembly_id.to_string()).or_default();
            current.push(item);
        }

        SequenceWithRegions {
            id: id.unwrap(),
            urs_taxid: urs_taxid.unwrap(),
            regions,
        }
    }
}

impl Iterator for RegionGrouper {
    type Item = SequenceWithRegions;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            None => None,
            Some(Grouped::Optional {
                id,
                data: _data,
            }) => {
                assert!(false, "Bad data at {}", id);
                None
            },
            Some(Grouped::Required {
                id,
                data: _data,
            }) => {
                assert!(false, "Bad data at {}", id);
                None
            },
            Some(Grouped::Multiple {
                id: _id,
                data,
            }) => Some(SequenceWithRegions::from_iter(data)),
        }
    }
}
