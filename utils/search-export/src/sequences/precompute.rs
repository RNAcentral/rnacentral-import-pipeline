use serde::{
    Deserialize,
    Serialize,
};
use std::{
    convert::TryFrom,
    path::Path,
};

use serde_with::CommaSeparator;

use anyhow::Result;
use phf::phf_map;
use rnc_core::grouper;

use crate::sequences::so_tree::SoId;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Precompute {
    pub id: usize,
    description: String,
    rna_type: String,
    has_coordinates: bool,
    so_rna_type: SoId,

    #[serde(with = "serde_with::rust::StringWithSeparator::<CommaSeparator>")]
    databases: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
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
    "SO:0001459" => "SO:0000655",
};

impl Precompute {
    pub fn so_rna_type(&self) -> SoId {
        // This is absurd
        let key: String = self.so_rna_type.clone().into();
        let k: &str = &key;
        let found = SO_MAPPING.get(k);
        if let Some(raw) = found {
            return SoId::try_from(raw.to_string()).unwrap();
        }
        self.so_rna_type.clone()
    }
}

impl From<&Precompute> for PrecomputeSummary {
    fn from(pre: &Precompute) -> PrecomputeSummary {
        Self {
            description: pre.description.to_owned(),
            rna_type: pre.rna_type.replace("_", " "),
            has_coordinates: pre.has_coordinates,
            databases: pre.databases.to_owned(),
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

impl PrecomputeSummary {
    /// Get a reference to the precompute summary's description.
    pub fn description(&self) -> &str {
        self.description.as_str()
    }

    pub fn databases(&self) -> &[String] {
        &self.databases
    }

    pub fn has_coordinates(&self) -> bool {
        self.has_coordinates
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn precompute(so_rna_type: &str) -> Precompute {
        serde_json::from_value(json!({
            "id": 1,
            "description": "a thing",
            "rna_type": "lnc_RNA",
            "has_coordinates": true,
            "so_rna_type": so_rna_type,
            "databases": "ENA,RefSeq",
        }))
        .unwrap()
    }

    fn so(id: &str) -> SoId {
        SoId::try_from(id).unwrap()
    }

    #[test]
    fn retired_so_terms_are_mapped_to_their_replacement() {
        assert_eq!(precompute("SO:0001263").so_rna_type(), so("SO:0000655"));
        assert_eq!(precompute("SO:0001637").so_rna_type(), so("SO:0000252"));
        assert_eq!(precompute("SO:0000122").so_rna_type(), so("SO:0000673"));
    }

    #[test]
    fn other_so_terms_are_left_alone() {
        assert_eq!(precompute("SO:0000374").so_rna_type(), so("SO:0000374"));
    }

    #[test]
    fn databases_are_split_on_commas() {
        let summary = PrecomputeSummary::from(&precompute("SO:0000374"));
        assert_eq!(summary.databases(), ["ENA", "RefSeq"]);
    }

    #[test]
    fn the_summary_shows_the_rna_type_with_spaces() {
        let summary = PrecomputeSummary::from(&precompute("SO:0000374"));
        assert_eq!(summary.rna_type, "lnc RNA");
        assert!(summary.has_coordinates());
        assert_eq!(summary.description(), "a thing");
    }

    #[test]
    fn coordinates_are_written_as_python_style_booleans() {
        assert_eq!(precompute("SO:0000374").has_coordiante_string(), "True");
    }
}
