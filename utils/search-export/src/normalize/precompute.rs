use serde::{
    Deserialize,
    Serialize,
};
use std::path::Path;

use serde_with::CommaSeparator;

use anyhow::Result;
use rnc_core::grouper;
use phf::phf_map;

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

static SO_MAPPING: phf::Map<&'static str, &'static str> = phf_map! {
    "SO:0001263" => "SO:0000655",
    "SO:0001637" => "SO:0000252",
    "SO:0001427" => "SO:0000673",
    "SO:1001268" => "SO:0000673",
    "SO:0000122" => "SO:0000673",
    "SO:0005836" => "SO:0000673",
    "SO:0000185" => "SO:0000673",
    "SO:0000204" => "SO:0000673",
    "SO:0000140" => "SO:0000673",
    "SO:0000243" => "SO:0000673",
    "SO:0000205" => "SO:0000673",
    "SO:1001274" => "SO:0000673",
    "SO:0000726" => "SO:0000673",
    "SO:0000077" => "SO:0000673",
    "SO:0000836" => "SO:0000673",
    "SO:0000233" => "SO:0000673",
    "SO:0001267" => "SO:0000275",
    "SO:0001243" => "SO:0001244",
};


impl Precompute {
    pub fn so_rna_type(&self) -> &str {
        let given: &str = &self.so_rna_type;
        SO_MAPPING.get(given).cloned().unwrap_or(given)
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

impl grouper::HasIndex for Precompute {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, max: usize, output: &Path) -> Result<()> {
    grouper::group::<Precompute>(grouper::Criteria::ExactlyOne, &path, 1, max, &output)
}
