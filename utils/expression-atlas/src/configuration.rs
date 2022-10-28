/// This module handles the parsing of the configuration file
use quick_xml::de::from_str;
use quick_xml::DeError;
use serde::Deserialize;
use std::fs;

use std::path::PathBuf;

#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Config {
    #[serde(rename = "experimentType")]
    pub exp_type: String,
    #[serde(rename = "analytics")]
    pub analytics: Vec<Analytics>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Analytics {
    pub assay_groups: AssayGroups,
    pub array_design: Option<String>,
    pub contrasts: Option<Contrasts>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct AssayGroups {
    pub assay_group: Vec<AssayGroup>,
}

#[derive(Debug, Deserialize, PartialEq, Eq, Clone)]
pub struct AssayGroup {
    pub id: String,
    pub label: Option<String>, // This contains the factors in a ; separated list
    #[serde(rename = "assay", default)]
    pub assays: Vec<String>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Contrasts {
    pub contrast: Vec<Contrast>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
pub struct Contrast {
    pub id: String,
    pub name: String,
    #[serde(alias = "reference_assay_group")]
    pub ref_group: String,
    #[serde(alias = "test_assay_group")]
    pub test_group: String,
}

pub fn parse_config(file: &PathBuf) -> Result<Config, DeError> {
    from_str::<Config>(&fs::read_to_string(file).unwrap())
}
