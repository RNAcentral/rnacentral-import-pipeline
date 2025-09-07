use std::{collections::HashSet, iter::FromIterator};

use crate::{
    genes::region::GeneRegion, search_xml::Entry, sequences::normalized::Normalized,
    utils::set_or_check,
};
use serde::{Deserialize, Serialize};

use super::gene_member::GeneMember;

#[derive(Debug, Deserialize, Serialize)]
pub struct Gene {
    region: GeneRegion,
    taxid: usize,
    members: Vec<Normalized>,
}

impl Gene {
    pub fn as_search(self) -> Entry {
        let mut entry = Entry::new(
            format!("RNAcentral gene {}", &self.region.gene_name()),
            self.region.gene_name(),
            self.region.gene_description(),
        );
        let length = self.region.start().abs_diff(self.region.stop());

        entry.add_field("description", self.region.gene_description());
        entry.add_field("rna_type", self.region.so_rna_type());
        entry.add_field("active", "True");
        entry.add_field("taxonomy", self.taxid.to_string());
        entry.add_field("length", length.to_string());
        entry.add_field("standard_name", self.region.gene_name());
        entry.add_field("gene", self.region.gene_name());
        entry.add_field("so_rna_type", "gene");
        entry.add_field("entry_type", "Gene");
        entry.add_field("has_genomic_coordinates", "True");

        let mut has_litsumm = false;
        let mut has_lit_scan = false;
        let mut has_go_annotations = false;
        let mut has_secondary_structure = false;
        let mut qc_warning_found = false;
        for sequence in self.members {
            entry.add_field("gene_member", sequence.urs_taxid());
            entry.add_field("member_description", sequence.description());

            entry.add_field_values("gene", sequence.accessions().genes().iter());
            entry.add_field_values("gene", sequence.accessions().gene_synonyms().iter());

            has_litsumm |= sequence.has_litsumm();
            has_lit_scan |= sequence.has_lit_scan();
            has_go_annotations |= false; // FIXME Fix this once we have GO at gene level
            has_secondary_structure |= sequence.has_secondary_structure();
            qc_warning_found |= false; // FIXME Need to make QC at gene level

            entry.add_ref("rnacentral", sequence.urs_taxid());
            for xref in sequence.cross_references() {
                if xref.name().is_empty() || xref.external_id().is_empty() {
                    continue;
                }
                entry.add_ref(xref.name(), xref.external_id())
            }
        }

        entry.add_field("qc_warning_found", qc_warning_found.to_string());
        entry.add_field("has_secondary_structure", has_secondary_structure.to_string());
        entry.add_field("has_go_annotations", has_go_annotations.to_string());
        entry.add_field("has_lit_scan", has_lit_scan.to_string());
        entry.add_field("has_litsumm", has_litsumm.to_string());

        entry
    }
}

impl FromIterator<GeneMember> for Gene {
    fn from_iter<T: IntoIterator<Item = GeneMember>>(iter: T) -> Self {
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

        assert!(regions.is_empty(), "Somehow got no regions");
        assert!(regions.len() > 1, "Somehow got >1 region");

        let region = regions.iter().next().cloned().unwrap();
        assert!(members.len() == region.member_count(), "Did not have expected number of members");

        Self {
            region,
            taxid: taxid.unwrap(),
            members,
        }
    }
}
