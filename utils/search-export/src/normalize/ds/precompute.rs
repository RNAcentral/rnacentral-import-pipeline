use serde::{
    Deserialize,
    Serialize,
};

use serde_with::CommaSeparator;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Precompute {
    pub id: usize,
    description: String,
    rna_type: String,
    has_coordinates: bool,
    so_rna_type: String,

    #[serde(with = "serde_with::rust::StringWithSeparator::<CommaSeparator>")]
    databases: Vec<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct PrecomputeSummary {
    description: String,
    rna_type: String,
    has_coordinates: bool,
    databases: Vec<String>,
}

impl Precompute {
    pub fn so_rna_type(&self) -> &str {
        &self.so_rna_type
    }
}

impl From<Precompute> for PrecomputeSummary {
    fn from(pre: Precompute) -> PrecomputeSummary {
        Self {
            description: pre.description,
            rna_type: pre.rna_type.replace("_", " "),
            has_coordinates: pre.has_coordinates,
            databases: pre.databases,
        }
    }
}

impl Precompute {
    pub fn has_coordiante_string(self) -> String {
        match self.has_coordinates {
            true => String::from("True"),
            false => String::from("False"),
        }
    }
}
