use serde::Serialize;

#[derive(Debug, Serialize, PartialEq)]
pub struct ValueOnly {
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct HierarchicalField {
    name: String,
    root: ValueOnly,
    #[serde(rename = "child")]
    children: Vec<ValueOnly>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct NamedField {
    name: String,

    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Serialize, PartialEq, Default)]
pub struct AdditionalFields {
    #[serde(rename = "field", default)]
    fields: Vec<NamedField>,

    #[serde(rename = "hierarchical_field", default)]
    trees: Vec<HierarchicalField>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct Ref {
    dbname: String,
    dbkey: String,
}

#[derive(Debug, Serialize, PartialEq, Default)]
pub struct CrossReferences {
    #[serde(rename = "ref", default)]
    refs: Vec<Ref>,
}

#[derive(Debug, Serialize, PartialEq)]
#[serde(rename = "entry")]
pub struct Entry {
    id: String,

    #[serde(rename = "name")]
    name: ValueOnly,

    #[serde(rename = "description")]
    description: ValueOnly,
    cross_references: CrossReferences,
    additional_fields: AdditionalFields,
}

impl<S> From<S> for ValueOnly
where
    S: Into<String>,
{
    fn from(value: S) -> Self {
        Self {
            value: value.into(),
        }
    }
}

// impl<'a> From<&'a str> for ValueOnly {
//     fn from(value: &'a str) -> Self {
//         Self {
//             value.to_string(),
//         }
//     }
// }
//
// impl<'a> From<&'a String> for ValueOnly {
//     fn from(value: &'a String) -> Self {
//         Self {
//             value,
//         }
//     }
// }
//
// impl<'a> ValueOnly<'a> {
//     pub fn new(value: &'a str) -> Self {
//         Self {
//             value,
//         }
//     }
// }

impl Entry {
    pub fn new<N, D>(id: String, name: N, description: D) -> Self
    where
        N: Into<String>,
        D: Into<String>,
    {
        Self {
            id,
            name: ValueOnly::from(name),
            description: ValueOnly::from(description),
            cross_references: CrossReferences::default(),
            additional_fields: AdditionalFields::default(),
        }
    }

    pub fn add_ref<D, K>(&mut self, database_name: D, key: K)
    where
        D: Into<String>,
        K: Into<String>,
    {
        self.cross_references.add_ref(database_name, key);
    }

    pub fn add_field<N, S>(&mut self, name: N, value: S)
    where
        N: Into<String>,
        S: Into<String>,
    {
        self.additional_fields.add_field(name, value);
    }

    pub fn add_field_values<N, C, V>(&mut self, name: N, values: V)
    where
        N: Into<String>,
        C: Into<String>,
        V: Iterator<Item = C>,
    {
        let name: String = name.into();
        for value in values {
            self.add_field(&name, value);
        }
    }

    pub fn add_tree<N, R, C, I>(&mut self, name: N, root: R, children: I)
    where
        N: Into<String>,
        R: Into<String>,
        C: Into<String>,
        I: Iterator<Item = C>,
    {
        self.additional_fields.add_tree(name, root, children);
    }
}

impl CrossReferences {
    pub fn add_ref<D, K>(&mut self, database_name: D, key: K)
    where
        D: Into<String>,
        K: Into<String>,
    {
        self.refs.push(Ref {
            dbname: database_name.into(),
            dbkey: key.into(),
        });
    }
}

impl AdditionalFields {
    pub fn add_field<N, V>(&mut self, name: N, value: V)
    where
        N: Into<String>,
        V: Into<String>,
    {
        self.fields.push(NamedField {
            name: name.into(),
            value: value.into(),
        });
    }

    pub fn add_tree<N, R, C, I>(&mut self, name: N, root: R, children: I)
    where
        N: Into<String>,
        R: Into<String>,
        C: Into<String>,
        I: Iterator<Item = C>,
    {
        let mut members = Vec::new();
        for child in children {
            members.push(ValueOnly::from(child));
        }
        self.trees.push(HierarchicalField {
            name: name.into(),
            root: root.into().into(),
            children: members,
        });
    }
}
