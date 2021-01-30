use std::path::Path;

use serde::{
    Deserialize,
    Serialize,
};

use itertools::Itertools;

use anyhow::Result;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use rnc_core::psql::JsonlIterator;

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub struct Accession {
    urs_taxid: String,
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
    urs_taxid: String,
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
    urs_taxid: String,
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
    urs_taxid: String,
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
    urs_taxid: String,
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
    urs_taxid: String,
    assembly_id: String,
    chromosome: String,
    strand: String,
    start: usize,
    stop: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Xref {
    urs_taxid: String,
    length: usize,
    deleted: String,
    last_release: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Metadata {
    upi: String,
    taxid: usize,
    length: usize,
    last_release: usize,
    coordinates: Vec<Coordinate>,
    deleted: bool,
    previous: Option<Previous>,
    rfam_hits: Vec<RfamHit>,
    r2dt_hits: Vec<R2dtHit>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Normalized {
    upi: String,
    taxid: usize,
    length: usize,
    last_release: usize,
    coordinates: Vec<Coordinate>,
    accessions: Vec<Accession>,
    deleted: bool,
    previous: Option<Previous>,
    rfam_hits: Vec<RfamHit>,
    r2dt_hits: Vec<R2dtHit>,
}

impl From<RawAccessionEntry> for Accession {
    fn from(raw: RawAccessionEntry) -> Self {
        return Self {
            urs_taxid: raw.urs_taxid,
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

impl Normalized {
    fn new(
        urs_taxid: String,
        raw_accessions: impl Iterator<Item = RawAccessionEntry>,
        raw_xrefs: impl Iterator<Item = Xref>,
        raw_coordinates: Option<impl Iterator<Item = Coordinate>>,
        raw_rfam_hits: Option<impl Iterator<Item = RfamHit>>,
        raw_r2dt_hits: Option<impl Iterator<Item = R2dtHit>>,
        raw_previous: Option<impl Iterator<Item = Previous>>,
    ) -> Result<Self> {
        let parts: Vec<&str> = urs_taxid.split("_").collect();
        let upi = parts[0].to_owned();
        let taxid = parts[1].parse::<usize>().unwrap();

        let xrefs: Vec<Xref> = raw_xrefs.collect();

        let coordinates = raw_coordinates.into_iter().flatten().collect();
        let rfam_hits = raw_rfam_hits.into_iter().flatten().collect();
        let r2dt_hits = raw_r2dt_hits.into_iter().flatten().collect();
        let previous: Vec<Previous> = raw_previous.into_iter().flatten().collect();
        let previous = match previous.len() {
            0 => None,
            _ => Some(previous[0].clone()),
        };

        let length = xrefs[0].length;
        let last_release = xrefs.iter().map(|x| x.last_release).max().unwrap();

        let mut deleted = true;
        let mut accessions = Vec::with_capacity(xrefs.len());
        for raw_accession in raw_accessions {
            if raw_accession.is_active {
                deleted = false;
            }
            accessions.push(Accession::from(raw_accession));
        }

        return Ok(Self {
            upi,
            taxid,
            length,
            last_release,
            coordinates,
            accessions,
            deleted,
            previous,
            rfam_hits,
            r2dt_hits,
        });
    }
}

impl Metadata {
    fn new(
        urs_taxid: String,
        raw_xrefs: impl Iterator<Item = Xref>,
        raw_coordinates: Option<impl Iterator<Item = Coordinate>>,
        raw_rfam_hits: Option<impl Iterator<Item = RfamHit>>,
        raw_r2dt_hits: Option<impl Iterator<Item = R2dtHit>>,
        raw_previous: Option<impl Iterator<Item = Previous>>,
    ) -> Result<Self> {
        let parts: Vec<&str> = urs_taxid.split("_").collect();
        let upi = parts[0].to_owned();
        let taxid = parts[1].parse::<usize>().unwrap();

        let xrefs: Vec<Xref> = raw_xrefs.collect();

        let coordinates = raw_coordinates.into_iter().flatten().collect();
        let rfam_hits = raw_rfam_hits.into_iter().flatten().collect();
        let r2dt_hits = raw_r2dt_hits.into_iter().flatten().collect();
        let previous: Vec<Previous> = raw_previous.into_iter().flatten().collect();
        let previous = match previous.len() {
            0 => None,
            _ => Some(previous[0].clone()),
        };

        let length = xrefs[0].length;
        let last_release = xrefs.iter().map(|x| x.last_release).max().unwrap();

        let deleted = xrefs.iter().all(|x| x.deleted == "Y");

        return Ok(Self {
            upi,
            taxid,
            length,
            last_release,
            coordinates,
            deleted,
            previous,
            rfam_hits,
            r2dt_hits,
        });
    }
}

pub fn write_metadata(
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    xref_file: &Path,
    output: &Path,
) -> Result<()> {
    let xrefs = JsonlIterator::from_path(xref_file)?;
    let xrefs = xrefs.group_by(|x: &Xref| x.urs_taxid.to_owned());
    let xrefs = xrefs.into_iter().assume_sorted_by_key();

    let coordinates = JsonlIterator::from_path(coordinate_file)?;
    let coordinates = coordinates.group_by(|c: &Coordinate| c.urs_taxid.to_owned());
    let coordinates = coordinates.into_iter().assume_sorted_by_key();

    let rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let rfam_hits = rfam_hits.group_by(|h: &RfamHit| h.urs_taxid.to_owned());
    let rfam_hits = rfam_hits.into_iter().assume_sorted_by_key();

    let r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let r2dt_hits = r2dt_hits.group_by(|h: &R2dtHit| h.urs_taxid.to_owned());
    let r2dt_hits = r2dt_hits.into_iter().assume_sorted_by_key();

    let previous = JsonlIterator::from_path(previous_file)?;
    let previous = previous.group_by(|p: &Previous| p.urs_taxid.to_owned());
    let previous = previous.into_iter().assume_sorted_by_key();

    let partial =
        xrefs.left_join(coordinates).left_join(rfam_hits).left_join(r2dt_hits).left_join(previous);

    let mut output = rnc_utils::buf_writer(output)?;
    for (id, data) in partial {
        let ((((xrefs, coordinates), rfam_hits), r2dt_hits), previous) = data;
        let norm = Metadata::new(id, xrefs, coordinates, rfam_hits, r2dt_hits, previous)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    Ok(())
}

pub fn write(
    accession_file: &Path,
    coordinate_file: &Path,
    rfam_hits_file: &Path,
    r2dt_hits_file: &Path,
    previous_file: &Path,
    xref_file: &Path,
    output: &Path,
) -> Result<()> {
    let accessions = JsonlIterator::from_path(accession_file)?;
    let accessions = accessions.group_by(|a: &RawAccessionEntry| a.urs_taxid.to_owned());
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let xrefs = JsonlIterator::from_path(xref_file)?;
    let xrefs = xrefs.group_by(|x: &Xref| x.urs_taxid.to_owned());
    let xrefs = xrefs.into_iter().assume_sorted_by_key();

    let coordinates = JsonlIterator::from_path(coordinate_file)?;
    let coordinates = coordinates.group_by(|c: &Coordinate| c.urs_taxid.to_owned());
    let coordinates = coordinates.into_iter().assume_sorted_by_key();

    let rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let rfam_hits = rfam_hits.group_by(|h: &RfamHit| h.urs_taxid.to_owned());
    let rfam_hits = rfam_hits.into_iter().assume_sorted_by_key();

    let r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let r2dt_hits = r2dt_hits.group_by(|h: &R2dtHit| h.urs_taxid.to_owned());
    let r2dt_hits = r2dt_hits.into_iter().assume_sorted_by_key();

    let previous = JsonlIterator::from_path(previous_file)?;
    let previous = previous.group_by(|p: &Previous| p.urs_taxid.to_owned());
    let previous = previous.into_iter().assume_sorted_by_key();

    let partial = accessions
        .join(xrefs)
        .left_join(coordinates)
        .left_join(rfam_hits)
        .left_join(r2dt_hits)
        .left_join(previous);

    let mut output = rnc_utils::buf_writer(output)?;
    for (id, data) in partial {
        let (((((accessions, xrefs), coordinates), rfam_hits), r2dt_hits), previous) = data;
        let norm =
            Normalized::new(id, accessions, xrefs, coordinates, rfam_hits, r2dt_hits, previous)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    Ok(())
}
