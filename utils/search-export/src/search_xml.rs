use serde::{
    Deserialize,
    Serialize,
};

use typed_builder::TypedBuilder;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Name {
    #[serde(rename = "$value")]
    value: String,
}

impl Name {
    pub fn new(value: String) -> Self { Self { value } }
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Description {
    #[serde(rename = "$value")]
    value: String,
}

impl Description {
    pub fn new(value: String) -> Self { Self { value } }
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct FieldChild {
    #[serde(rename = "$value")]
    value: String
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct HierarchicalField {
    name: String,

    root: FieldChild,
    #[serde(rename = "child")]
    children: Vec<FieldChild>
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
            trees: Vec::new()
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
    pub fn add_field(&mut self, name: &str, value: &str) {
        self.fields.push(ValueField {
            name: name.to_string(),
            value: value.to_string(),
        });
    }

    pub fn add_tree(&mut self, name: &str, root: FieldChild, children: Vec<FieldChild>) {
        self.trees.push(HierarchicalField {
            name: name.to_string(),
            root,
            children,
        });
    }
}
