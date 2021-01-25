use std::{
    convert::From,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Accession {
    accession: String,
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
    all_species: Vec<String>,
    all_common_names: Vec<String>,
    so_rna_type: Option<String>,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct RfamHit {
    rfam_hit_id: usize,
    model: String,
    model_rna_type: Option<String>,
    model_domain: Option<String>,
    model_name: String,
    model_long_name: String,
    model_completeness: f64,
    model_start: usize,
    model_stop: usize,
    sequence_completeness: f64,
    sequence_start: usize,
    sequence_stop: usize,
}

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct R2dtHit {
    model_id: String,
    model_name: String,
    model_source: String,
    model_so_term: Option<String>,
    sequence_coverage: Option<f64>,
    model_coverage: Option<f64>,
    sequence_basepairs: Option<usize>,
    model_basepairs: Option<usize>,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct Previous {
    id: String,
    upi: String,
    taxid: usize,
    databases: String,
    description: Option<String>,
    has_coordinates: bool,
    is_active: bool,
    is_fragment: bool,
    is_locus_representative: bool,
    last_release: usize,
    rfam_problems: String,
    rna_type: String,
    short_description: Option<String>,
    so_rna_type: Option<String>,
    update_date: Option<String>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct RawAccessionEntry {
    accession: String,
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

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Coordinate {
    assembly_id: String,
    chromosome: String,
    strand: String,
    start: usize,
    stop: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Xref {
    length: usize,
    deleted: String,
    last_release: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Raw {
    id: String,
    #[serde(rename = "xref")]
    xrefs: Vec<Xref>,
    accessions: Vec<RawAccessionEntry>,
    coordinates: Vec<Coordinate>,

    #[serde(rename = "rfam-hits")]
    rfam_hits: Vec<RfamHit>,

    #[serde(rename = "r2dt-hits")]
    r2dt_hits: Vec<R2dtHit>,
    previous: Vec<Previous>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Normalized {
    upi: String,
    taxid: usize,
    length: usize,
    last_release: usize,
    coordinates: Vec<Coordinate>,
    accessions: Vec<Accession>,
    deleted: Vec<bool>,
    previous: Option<Previous>,
    rfam_hits: Vec<RfamHit>,
    r2dt_hits: Vec<R2dtHit>,
}

impl From<RawAccessionEntry> for Accession {
    fn from(raw: RawAccessionEntry) -> Self {
        return Self {
            accession: raw.accession,
            is_active: raw.is_active,
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

impl From<Raw> for Normalized {
    fn from(raw: Raw) -> Self {
        let parts: Vec<&str> = raw.id.split("_").collect();
        let upi = parts[0].to_owned();
        let taxid = parts[1].parse::<usize>().unwrap();
        let previous = raw.previous.get(0).cloned();
        let length = raw.xrefs[0].length;
        let last_release = raw.xrefs.iter().map(|x| x.last_release).max().unwrap();

        let mut deleted = true;
        let mut accessions = Vec::with_capacity(raw.accessions.len());
        for raw_accession in raw.accessions {
            if raw_accession.is_active {
                deleted = false;
            }
            accessions.push(Accession::from(raw_accession));
        }

        return Self {
            upi,
            taxid,
            length,
            last_release,
            coordinates: raw.coordinates,
            accessions,
            deleted: vec![deleted],
            previous,
            rfam_hits: raw.rfam_hits,
            r2dt_hits: raw.r2dt_hits,
        };
    }
}

pub fn write(filename: &Path, output: &Path) -> Result<()> {
    let mut reader = rnc_utils::buf_reader(filename)?;
    let mut output = rnc_utils::buf_writer(output)?;
    let mut buf = String::new();
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                let cleaned = buf.replace("\\\\", "\\");
                let raw: Raw = serde_json::from_str(&cleaned)?;
                let norm: Normalized = Normalized::from(raw);
                serde_json::to_writer(&mut output, &norm)?;
                writeln!(&mut output)?;
                buf.clear();
            },
        }
    }

    Ok(())
}
