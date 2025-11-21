use serde::{
    Deserialize,
    Serialize,
};

use std::path::Path;

use anyhow::Result;

use rnc_core::grouper;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Accession {
    pub urs_taxid: String,
    pub accession: String,
    pub is_active: bool,
    pub last_release: usize,
    description: String,
    gene: Option<String>,
    optional_id: Option<String>,
    database: String,
    species: Option<String>,
    common_name: Option<String>,
    feature_name: Option<String>,
    ncrna_class: Option<String>,
    locus_tag: Option<String>,
    organelle: Option<String>,
    lineage: Option<String>,
    all_species: Vec<String>,
    all_common_names: Vec<String>,
    so_rna_type: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RawAccessionEntry {
    pub id: usize,
    pub urs_id: usize,
    pub urs_taxid: String,
    pub accession: String,
    last_release: usize,
    is_active: bool,
    description: String,
    gene: Option<String>,
    optional_id: Option<String>,
    database: String,
    species: Option<String>,
    common_name: Option<String>,
    feature_name: Option<String>,
    ncrna_class: Option<String>,
    locus_tag: Option<String>,
    organelle: Option<String>,
    lineage: Option<String>,
    all_species: Vec<Option<String>>,
    all_common_names: Vec<Option<String>>,
    so_rna_type: Option<String>,
}

impl From<RawAccessionEntry> for Accession {
    fn from(raw: RawAccessionEntry) -> Self {
        return Self {
            urs_taxid: raw.urs_taxid,
            accession: raw.accession,
            is_active: raw.is_active,
            last_release: raw.last_release,
            description: raw.description,
            gene: raw.gene,
            optional_id: raw.optional_id,
            database: raw.database,
            species: raw.species,
            common_name: raw.common_name,
            feature_name: raw.feature_name,
            ncrna_class: raw.ncrna_class,
            locus_tag: raw.locus_tag,
            organelle: raw.organelle,
            lineage: raw.lineage,
            all_species: raw.all_species.into_iter().filter_map(|s| s).collect(),
            all_common_names: raw.all_common_names.into_iter().filter_map(|s| s).collect(),
            so_rna_type: raw.so_rna_type,
        };
    }
}

impl grouper::HasIndex for RawAccessionEntry {
    fn index(&self) -> usize {
        self.id
    }
}

pub fn group(path: &Path, min: usize, max: usize, output: &Path) -> Result<()> {
    grouper::group::<RawAccessionEntry>(grouper::Criteria::AtleastOne, &path, min, max, &output)
}
