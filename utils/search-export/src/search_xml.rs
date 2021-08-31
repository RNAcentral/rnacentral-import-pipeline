use serde::Serialize;

use typed_builder::TypedBuilder;

#[derive(Debug, Serialize, PartialEq)]
pub struct ValueOnly<'a> {
    #[serde(rename = "$value")]
    value: &'a str,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct HierarchicalField<'a> {
    name: &'a str,
    root: ValueOnly<'a>,
    #[serde(rename = "child")]
    children: Vec<ValueOnly<'a>>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct NamedField<'a> {
    name: &'a str,
    #[serde(rename = "$value")]
    value: &'a str,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct AdditionalFields<'a> {
    #[serde(rename = "field", default)]
    fields: Vec<NamedField<'a>>,

    #[serde(rename = "hierarchical_field", default)]
    trees: Vec<HierarchicalField<'a>>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct Ref<'a> {
    dbname: &'a str,
    dbkey: &'a str,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct CrossReferences<'a> {
    #[serde(rename = "ref", default)]
    refs: Vec<Ref<'a>>,
}

#[derive(TypedBuilder, Debug, Serialize, PartialEq)]
pub struct Entry<'a> {
    id: &'a str,
    name: ValueOnly<'a>,
    description: ValueOnly<'a>,
    cross_references: CrossReferences<'a>,
    additional_fields: AdditionalFields<'a>,
}

impl<'a> Default for CrossReferences<'a> {
    fn default() -> Self {
        Self {
            refs: Vec::new(),
        }
    }
}

impl<'a> Default for AdditionalFields<'a> {
    fn default() -> Self {
        Self {
            fields: Vec::new(),
            trees: Vec::new(),
        }
    }
}

impl<'a> From<&'a str> for ValueOnly<'a> {
    fn from(value: &'a str) -> Self {
        Self {
            value,
        }
    }
}

impl<'a> From<&'a String> for ValueOnly<'a> {
    fn from(value: &'a String) -> Self {
        Self {
            value,
        }
    }
}

impl<'a> ValueOnly<'a> {
    pub fn new(value: &'a str) -> Self {
        Self {
            value,
        }
    }
}

impl<'a> Entry<'a> {
    pub fn new(id: &'a str, name: &'a str, description: &'a str) -> Self {
        Self {
            id,
            name: ValueOnly::new(name),
            description: ValueOnly::new(description),
            cross_references: CrossReferences::default(),
            additional_fields: AdditionalFields::default(),
        }
    }

    pub fn add_ref(&mut self, database_name: &'a str, key: &'a str) {
        self.cross_references.add_ref(database_name, key);
    }

    pub fn add_field(&mut self, name: &'a str, value: &'a str) {
        self.additional_fields.add_field(name, value);
    }

    pub fn add_tree(&mut self, name: &'a str, root: &'a str, children: &'a Vec<String>) {
        self.additional_fields.add_tree(name, root, children);
    }
}

impl<'a> CrossReferences<'a> {
    pub fn add_ref(&mut self, database_name: &'a str, key: &'a str) {
        self.refs.push(Ref {
            dbname: database_name,
            dbkey: key,
        });
    }
}

impl<'a> AdditionalFields<'a> {
    pub fn add_field(&mut self, name: &'a str, value: &'a str) {
        self.fields.push(NamedField {
            name,
            value,
        });
    }

    pub fn add_tree(&mut self, name: &'a str, root: &'a str, children: &'a Vec<String>) {
        self.trees.push(HierarchicalField {
            name,
            root: root.into(),
            children: children.iter().map(|s| s.into()).collect(),
        });
    }
}
