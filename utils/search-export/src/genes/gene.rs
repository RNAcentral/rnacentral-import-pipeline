use std::iter::FromIterator;

use crate::{
    search_xml::{
        AdditionalFields,
        CrossReferences,
        Description,
        Entry,
        Name,
    },
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

impl Into<Entry> for Gene {
    fn into(self) -> Entry {
        let mut cross_references = CrossReferences::default();
        let mut additional_fields = AdditionalFields::default();

        for member in self.members {
            additional_fields.add_field("entry_type", "Gene");
            additional_fields.add_field("gene_member", member.urs_taxid());
            additional_fields.add_field("member_description", member.description());

            cross_references.add_ref("rnacentral", member.urs_taxid());
            for xref in member.cross_references() {
                cross_references.add_ref(xref.name(), xref.external_id())
            }
        }

        Entry::builder()
            .id(self.name.to_string())
            .name(Name::new("".to_string()))
            .description(Description::new("".to_string()))
            .cross_references(cross_references)
            .additional_fields(additional_fields)
            .build()
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
