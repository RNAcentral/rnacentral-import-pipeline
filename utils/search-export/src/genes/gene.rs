use std::iter::FromIterator;

use crate::{search_xml::Entry, utils::set_or_check};
use serde::{Deserialize, Serialize};

use super::gene_member::GeneMember;

#[derive(Debug, Deserialize, Serialize)]
pub struct Gene {
    id: usize,
    name: String,
    assembly_id: String,
    member_count: usize,
    members: Vec<GeneMember>,
}

impl<'a> Gene {
    pub fn fill_default_name(&mut self) {
        if self.name.is_empty() {
            self.name = format!("RNA Gene {}", &self.id)
        }
    }

    pub fn as_search(&'a self) -> Entry<'a> {
        let description = self.members[0].description();
        let mut entry = Entry::new(self.id.to_string(), &self.name, &description);

        entry.add_field("description", "TODO");
        entry.add_field("rna_type", "TODO");
        entry.add_field("active", "True");
        entry.add_field("gene", "TODO");
        entry.add_field("taxonomy", "TODO");
        entry.add_field("lineage_path", "TODO");
        entry.add_field("length", "TODO");
        entry.add_tree("so_rna_type", "TODO", &vec!["TODO".to_string()]);

        let mut has_litsumm = false;
        let mut has_lit_scan = false;
        let mut has_go_annotations = false;
        let mut has_secondary_structure = false;
        let mut qc_warning_found = false;
        let mut has_genomic_coordinates = false;
        for member in &self.members {
            entry.add_field("entry_type", "Gene");
            entry.add_field("gene_member", member.urs_taxid());
            entry.add_field("member_description", member.description());

            has_litsumm = has_litsumm | member.has_litsumm();
            has_lit_scan = has_lit_scan | member.has_lit_scan();
            has_go_annotations = has_go_annotations | member.has_go_annotations();
            has_secondary_structure = has_secondary_structure | member.has_secondary_structure();
            qc_warning_found = qc_warning_found | member.qc_warning_found();
            has_genomic_coordinates = has_genomic_coordinates | member.has_genomic_coordinates();

            entry.add_ref("rnacentral", member.urs_taxid());
            for xref in member.cross_references() {
                if xref.name().is_empty() || xref.external_id().is_empty() {
                    continue;
                }
                entry.add_ref(xref.name(), xref.external_id())
            }
        }

        entry.add_field("has_genomic_coordinates", &has_litsumm.to_string());
        entry.add_field("qc_warning_found", &qc_warning_found.to_string());
        entry.add_field("has_secondary_structure", &has_secondary_structure.to_string());
        entry.add_field("has_go_annotations", &has_go_annotations.to_string());
        entry.add_field("has_lit_scan", &has_lit_scan.to_string());
        entry.add_field("has_litsumm", &has_litsumm.to_string());

        entry
    }
}

impl FromIterator<GeneMember> for Gene {
    fn from_iter<T: IntoIterator<Item = GeneMember>>(iter: T) -> Self {
        let mut locus_id = None;
        let mut name = None;
        let mut member_count = None;
        let mut assembly_id = None;
        let mut members = Vec::new();

        for locus in iter {
            set_or_check(&mut locus_id, locus.locus_id());
            set_or_check(&mut name, locus.name().to_string());
            set_or_check(&mut member_count, *locus.member_count());
            set_or_check(&mut assembly_id, locus.assembly_id().to_string());
            members.push(GeneMember::from(locus));
        }

        Self {
            id: locus_id.unwrap(),
            name: name.unwrap().to_string(),
            member_count: member_count.unwrap().to_owned(),
            assembly_id: assembly_id.unwrap().to_string(),
            members,
        }
    }
}
