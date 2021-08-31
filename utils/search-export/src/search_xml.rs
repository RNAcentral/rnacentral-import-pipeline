use serde::{
    Deserialize,
    Serialize,
};

use typed_builder::TypedBuilder;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct ValueOnly {
    #[serde(rename = "$value")]
    value: String
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct HierarchicalField {
    name: String,

    root: ValueOnly,
    #[serde(rename = "child")]
    children: Vec<ValueOnly>
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct ValueField {
    name: String,
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct AdditionalFields {
    #[serde(rename = "field", default)]
    fields: Vec<ValueField>,

    #[serde(rename = "hierarchical_field", default)]
    trees: Vec<HierarchicalField>
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

#[derive(TypedBuilder, Debug, Deserialize, Serialize, PartialEq)]
pub struct Entry {
    id: String,
    name: ValueOnly,
    description: ValueOnly,
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
            trees: Vec::new()
        }
    }
}

impl From<String> for ValueOnly {
    fn from(value: String) -> Self {
        Self { value }
    }
}

impl ValueOnly {
    pub fn new(value: String) -> Self { Self { value } }
}

impl Entry {
    pub fn new(id: String, name: String, description: String) -> Self {
        Self {
            id,
            name: ValueOnly::new(name),
            description: ValueOnly::new(description),
            cross_references: CrossReferences::default(),
            additional_fields: AdditionalFields::default(),
        }
    }

    pub fn add_ref(&mut self, database_name: &str, key: &str) {
        self.cross_references.add_ref(database_name, key);
    }

    pub fn add_field(&mut self, name: &str, value: &str) {
        self.additional_fields.add_field(name, value);
    }

    pub fn add_tree(&mut self, name: &str, root: &str, children: Vec<String>) {
        self.additional_fields.add_tree(name, root, children);
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
    pub fn add_field(&mut self, name: &str, value: &str) {
        self.fields.push(ValueField {
            name: name.to_string(),
            value: value.to_string(),
        });
    }

    pub fn add_tree(&mut self, name: &str, root: &str, children: Vec<String>) {
        self.trees.push(HierarchicalField {
            name: name.to_string(),
            root: root.to_string().into(),
            children: children.iter().map(|s| s.to_string().into()).collect(),
        });
    }
}
