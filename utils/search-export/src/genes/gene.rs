use std::iter::FromIterator;

use crate::{
    search_xml::Entry,
    utils::set_or_check,
};
use serde::{
    Deserialize,
    Serialize,
};

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
        let mut entry = Entry::new(self.id.to_string(), &self.name, description);

        for member in &self.members {
            entry.add_field("entry_type", "Gene");
            entry.add_field("gene_member", member.urs_taxid());
            entry.add_field("member_description", member.description());

            entry.add_ref("rnacentral", member.urs_taxid());
            for xref in member.cross_references() {
                if xref.name().is_empty() || xref.external_id().is_empty() {
                    continue;
                }
                entry.add_ref(xref.name(), xref.external_id())
            }
        }

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
