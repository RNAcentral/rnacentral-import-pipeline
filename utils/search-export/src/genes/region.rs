use std::{
    collections::HashMap,
    iter::FromIterator,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;
use rnc_core::grouper::{
    self,
    Grouped,
};

use crate::{
    sequences::so_tree::SoId,
    utils::set_or_check,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct SequenceWithRegions {
    id: usize,
    urs_taxid: String,
    regions: HashMap<String, Vec<UrsRegion>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Deserialize, Serialize)]
pub struct GeneRegion {
    gene_id: usize,
    gene_name: String,
    gene_description: String,
    so_rna_type: SoId,
    short_description: String,
    start: usize,
    stop: usize,
    assembly_id: String,
    member_count: usize,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct UrsRegion {
    id: usize,
    urs_taxid: String,
    gene_id: usize,
    gene_name: String,
    gene_description: String,
    so_rna_type: SoId,
    short_description: String,
    start: usize,
    stop: usize,
    assembly_id: String,
    member_count: usize,
}

impl From<UrsRegion> for GeneRegion {
    fn from(value: UrsRegion) -> Self {
        Self {
            gene_id: value.gene_id,
            gene_name: value.gene_name,
            gene_description: value.gene_description,
            so_rna_type: value.so_rna_type,
            short_description: value.short_description,
            start: value.start,
            stop: value.stop,
            assembly_id: value.assembly_id,
            member_count: value.member_count,
        }
    }
}

impl UrsRegion {
    /// Get a reference to the region's id.
    pub fn id(&self) -> &usize {
        &self.id
    }

    /// Get a reference to the region's urs taxid.
    pub fn urs_taxid(&self) -> &str {
        self.urs_taxid.as_str()
    }

    pub fn gene_id(&self) -> usize {
        self.gene_id
    }
}

impl GeneRegion {
    pub fn gene_id(&self) -> usize {
        self.gene_id
    }

    pub fn gene_name(&self) -> &str {
        &self.gene_name
    }

    pub fn gene_description(&self) -> &str {
        &self.gene_description
    }

    pub fn so_rna_type(&self) -> &SoId {
        &self.so_rna_type
    }

    pub fn short_description(&self) -> &str {
        &self.short_description
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn stop(&self) -> usize {
        self.stop
    }

    pub fn assembly_id(&self) -> &str {
        &self.assembly_id
    }

    pub fn member_count(&self) -> usize {
        self.member_count
    }
}

pub struct RegionGrouper<T>
where
    T: Iterator<Item = Grouped<UrsRegion>>,
{
    iter: T,
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<UrsRegion>(grouper::Criteria::AnyNumber, path, 1, max, output)
}

impl grouper::HasIndex for UrsRegion {
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
    pub fn regions_by_assembly(&self) -> &HashMap<String, Vec<UrsRegion>> {
        &self.regions
    }
}

impl<T> RegionGrouper<T>
where
    T: Iterator<Item = Grouped<UrsRegion>>,
{
    pub fn new(iter: T) -> Self {
        Self {
            iter,
        }
    }
}

impl FromIterator<UrsRegion> for SequenceWithRegions {
    fn from_iter<T: IntoIterator<Item = UrsRegion>>(iter: T) -> Self {
        let mut regions: HashMap<String, Vec<UrsRegion>> = HashMap::new();
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

impl<T> Iterator for RegionGrouper<T>
where
    T: Iterator<Item = Grouped<UrsRegion>>,
{
    type Item = SequenceWithRegions;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iter.next() {
            None => None,
            Some(Grouped::Optional {
                id,
                data: _data,
            }) => {
                panic!("Bad data at {}", id);
            },
            Some(Grouped::Required {
                id,
                data: _data,
            }) => {
                panic!("Bad data at {}", id);
            },
            Some(Grouped::Multiple {
                id: _id,
                data,
            }) => Some(SequenceWithRegions::from_iter(data)),
        }
    }
}
