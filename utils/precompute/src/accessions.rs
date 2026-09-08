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

/// AnyNumber, not AtleastOne: a pair whose xrefs are gone has no accessions,
/// and failing here would abandon the other 25,000 pairs in the range.
/// Normalize drops the empty groups.
pub fn group(path: &Path, min: usize, max: usize, output: &Path) -> Result<()> {
    grouper::group::<RawAccessionEntry>(grouper::Criteria::AnyNumber, &path, min, max, &output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        path::PathBuf,
        sync::atomic::{
            AtomicUsize,
            Ordering as AtomicOrdering,
        },
    };

    fn temp_path(name: &str) -> PathBuf {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
        std::env::temp_dir().join(format!(
            "precompute-accessions-{}-{}-{}",
            std::process::id(),
            n,
            name
        ))
    }

    fn raw(id: usize) -> String {
        format!(
            r#"{{"id":{id},"urs_id":{id},"urs_taxid":"URS000000000{id}_9606","accession":"A{id}","last_release":1,"is_active":true,"description":"An RNA","gene":null,"optional_id":null,"database":"ENA","species":null,"common_name":null,"feature_name":null,"ncrna_class":null,"locus_tag":null,"organelle":null,"lineage":null,"all_species":[],"all_common_names":[],"so_rna_type":null}}"#,
            id = id
        )
    }

    /// A urs_taxid with no accessions at all (its xrefs are gone, so
    /// insert-chunk.sql built nothing) must group to an empty entry rather than
    /// failing the whole range with "Missing data for id".
    #[test]
    fn allows_ids_without_any_accessions() -> Result<()> {
        let input = temp_path("raw.json");
        let output = temp_path("grouped.json");
        fs::write(&input, format!("{}\n{}\n", raw(1), raw(3)))?;

        group(&input, 1, 4, &output)?;

        let written = fs::read_to_string(&output)?;
        let lines: Vec<&str> = written.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[1], r#"{"Multiple":{"id":2,"data":[]}}"#);

        fs::remove_file(&input)?;
        fs::remove_file(&output)?;
        Ok(())
    }
}
