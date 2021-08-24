use serde::{
    Deserialize,
    Serialize,
};

use crate::genes::{
    info::GeneInfo,
    member::GeneMembers,
};

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Name {
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Description {
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Field {
    name: String,
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct AdditionalFields {
    #[serde(rename = "field", default)]
    fields: Vec<Field>,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Ref {
    dbname: String,
    dbkey: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct CrossReferences {
    #[serde(rename = "ref", default)]
    refs: Vec<Ref>,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Entry {
    id: String,
    name: Name,
    description: Description,
    cross_references: CrossReferences,
    additional_fields: AdditionalFields,
}

impl Default for CrossReferences {
    fn default() -> Self {
        Self {
            refs: Vec::new(),
        }
    }
}

impl Default for AdditionalFields {
    fn default() -> Self {
        Self {
            fields: Vec::new(),
        }
    }
}

impl CrossReferences {
    pub fn add_ref(&mut self, database_name: &str, key: &str) {
        self.refs.push(Ref {
            dbname: database_name.to_string(),
            dbkey: key.to_string(),
        });
    }
}

impl AdditionalFields {
    pub fn add_field(&mut self, name: &str, value: String) {
        self.fields.push(Field {
            name: name.to_string(),
            value,
        });
    }
}

impl Entry {
    pub fn new(gene: GeneInfo, members: GeneMembers) -> Self {
        let mut cross_references = CrossReferences::default();
        let mut additional_fields = AdditionalFields::default();

        for sequence in members.sequences() {
            additional_fields.add_field("entry_type", "Gene".to_string());
            additional_fields.add_field("gene_member", sequence.urs_taxid().to_string());
            cross_references.add_ref("rnacentral", sequence.urs_taxid());
            for xref in sequence.cross_references() {
                cross_references.add_ref(xref.name(), xref.external_id())
            }
        }

        Self {
            id: gene.name().to_string(),
            name: Name {
                value: "".to_string(),
            },
            description: Description {
                value: "".to_string(),
            },
            cross_references,
            additional_fields,
        }
    }
}
