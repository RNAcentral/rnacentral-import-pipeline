use std::{collections::HashSet, iter::FromIterator};

use crate::{
    fields::{GeneEntry, GeneFields, SoRnaTreeField},
    genes::region::GeneRegion,
    search_xml::{SearchEntry, SearchValue},
    sequences::{
        accession::CrossReference,
        normalized::Normalized,
        so_tree::{SoMapping, SoName, SoTree},
    },
    utils::set_or_check,
};
use serde::{Deserialize, Serialize};

use super::gene_member::GeneMember;
use log;

#[derive(Debug, Deserialize, Serialize)]
pub struct Gene {
    region: GeneRegion,
    taxid: usize,
    so_tree: SoTree,
    members: Vec<Normalized>,
}

impl Gene {
    pub fn new<T: IntoIterator<Item = GeneMember>>(
        iter: T,
        trees: &SoMapping,
    ) -> anyhow::Result<Self> {
        let mut taxid = None;
        let mut members = Vec::new();
        let mut regions: HashSet<GeneRegion> = HashSet::new();

        for locus in iter {
            let (region, sequence) = locus.into_inner();
            let (_, tid) = sequence.urs_taxid().split_once('_').expect("Invalid URS_taxid");
            set_or_check(&mut taxid, tid.parse::<usize>().expect("Invalid taxid"));
            regions.insert(region.into());
            members.push(sequence);
        }

        assert!(!members.is_empty(), "Did not have any members");
        assert!(!regions.is_empty(), "Somehow got no regions for: {:?}", members.first());
        assert!(regions.len() == 1, "Somehow got >1 region for: {:?}", members.first());

        // Get the first region, this is safe at this point.
        let region = regions.iter().next().cloned().unwrap();

        if members.len() != region.member_count() {
            log::warn!(
                "Expected to find {} members, but found: {}",
                region.member_count(),
                members.len()
            );
        }

        let rna_type = region.so_rna_type();
        let so_tree = trees
            .tree(region.so_rna_type())
            .ok_or_else(|| anyhow::anyhow!("Failed to get SO tree for {rna_type}"))?;

        Ok(Self {
            region,
            taxid: taxid.unwrap(),
            so_tree,
            members,
        })
    }

    pub fn expert_databases(&self) -> HashSet<&String> {
        let mut expert_dbs = HashSet::new();
        for member in &self.members {
            expert_dbs.extend(member.precompute_summary().databases());
        }
        expert_dbs
    }

    pub fn length(&self) -> usize {
        self.region.start().abs_diff(self.region.stop())
    }

    pub fn so_type_tree(&self) -> (SoName, Vec<SoName>) {
        let mut names: Vec<SoName> = self.so_tree.names().into_iter().cloned().collect();
        let root = names.remove(0);
        (root, names)
    }

    pub fn so_name(&self) -> SoName {
        self.so_tree.term_name().unwrap().clone()
    }

    pub fn has_litsumm(&self) -> bool {
        self.members.iter().any(|m| m.has_litsumm())
    }

    pub fn has_lit_scan(&self) -> bool {
        self.members.iter().any(|m| m.has_lit_scan())
    }

    pub fn boost(&self) -> usize {
        let mut boost = 0usize;
        if self.taxid == 9606 {
            boost += 1;
        }
        boost
    }

    pub fn member_urs_taxids(&self) -> impl Iterator<Item = &str> {
        self.members.iter().map(|m| m.urs_taxid())
    }
}

impl SearchEntry<GeneEntry> for Gene {
    fn id(&self) -> &str {
        self.region.gene_name()
    }

    fn name(&self) -> &str {
        self.region.gene_name()
    }

    fn description(&self) -> &str {
        self.region.gene_description()
    }

    fn taxid(&self) -> usize {
        self.taxid
    }

    fn field_values(&self, field: &GeneFields) -> SearchValue {
        match field {
            GeneFields::Length => self.length().into(),
            GeneFields::ExpertDb => self.expert_databases().into(),
            GeneFields::Gene => self.region.gene_name().into(),
            GeneFields::GeneSynonym => SearchValue::empty(),
            GeneFields::HasGoAnnotations => false.into(),
            GeneFields::HasGenomicCoordinates => true.into(),
            GeneFields::Boost => self.boost().into(),
            GeneFields::StandardName => self.region.gene_name().into(),
            GeneFields::QcWarningFound => false.into(),
            GeneFields::GeneMember => SearchValue::from_iter(self.member_urs_taxids()),
            GeneFields::HasSecondaryStructure => false.into(),
            GeneFields::HasLitScan => self.has_lit_scan().into(),
            GeneFields::HasLitsumm => self.has_litsumm().into(),
            GeneFields::HasEditingEvent => false.into(),
            GeneFields::InsdcRNAType => "Gene".into(),
            GeneFields::SoRnaTypeName => self.so_name().into(),
            GeneFields::PubliGeneName => self.region.gene_name().into(),
        }
    }

    fn cross_references(&self) -> impl IntoIterator<Item = CrossReference> {
        let mut xrefs = Vec::new();
        for member in &self.members {
            xrefs.extend(member.cross_references().clone());
        }
        xrefs
    }

    fn tree_field_values(&self, field: &SoRnaTreeField) -> (String, Vec<String>) {
        match field {
            SoRnaTreeField::SoRnaType => {
                let (root, children) = self.so_type_tree();
                (root.to_string(), children.into_iter().map(|c| c.to_string()).collect())
            },
        }
    }
}
