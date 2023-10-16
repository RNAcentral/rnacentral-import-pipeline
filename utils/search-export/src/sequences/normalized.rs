use std::{
    convert::TryInto,
    iter::FromIterator,
    num::TryFromIntError,
};

use serde::{
    Deserialize,
    Serialize,
};

use thiserror::Error;

use rnc_core::{
    urs,
    urs_taxid,
};

use crate::sequences::{
    accession::{
        AccessionVec,
        CrossReference,
        RawAccession,
        ReferenceVec,
    },
    basic::Basic,
    crs::CrsVec,
    editing_events::EditingEvent,
    feedback::FeedbackVec,
    go_annotation::GoAnnotation,
    interacting_protein::InteractingProtein,
    interacting_rna::InteractingRna,
    litsumm::LitsummSummaries,
    orf::OrfVec,
    precompute::PrecomputeSummary,
    qa_status::QaStatus,
    r2dt::R2dt,
    raw::Raw,
    rfam_hit::RfamHitVec,
    so_tree,
};

#[derive(Error, Debug)]
pub enum NormalizationError {
    #[error("Could not parse {0}")]
    UrsParsingError(#[from] urs::Error),

    #[error("Could not parse {0}")]
    UrsTaxidParsingError(#[from] urs_taxid::Error),

    #[error("Could not convert {0} to urs number")]
    UrsNumberError(#[from] TryFromIntError),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Normalized {
    urs: String,
    taxid: usize,
    short_urs: String,
    deleted: String,
    qa_status: QaStatus,
    secondary: Option<R2dt>,
    cross_references: Vec<CrossReference>,
    crs: CrsVec,
    overlaps: FeedbackVec,
    go_annotations: Vec<GoAnnotation>,
    interacting_proteins: Vec<InteractingProtein>,
    interacting_rnas: Vec<InteractingRna>,
    publication_count: usize,
    litsumm: Vec<LitsummSummaries>,
    editing_events: Vec<EditingEvent>,
    so_rna_type_tree: so_tree::SoTree,

    #[serde(flatten)]
    orfs: OrfVec,

    #[serde(flatten)]
    pre_summary: PrecomputeSummary,

    #[serde(flatten)]
    basic: Basic,

    #[serde(flatten)]
    accessions: AccessionVec,

    #[serde(flatten)]
    references: ReferenceVec,

    #[serde(flatten)]
    rfam_hits: RfamHitVec,
    /* #[serde(flatten)]
     * dates: Dates, */
}

impl Normalized {
    pub fn new(raw: Raw, accessions: Vec<RawAccession>) -> Result<Self, NormalizationError> {
        let base = raw.base().clone();

        let urs_taxid = base.urs_taxid.to_owned();
        assert!(
            urs_taxid == accessions[0].urs_taxid,
            "Accession urs {} and base urs {} disagree for id {}",
            accessions[0].urs_taxid,
            urs_taxid,
            base.id
        );
        let parsed: urs_taxid::UrsTaxid = urs_taxid.parse()?;
        let references = ReferenceVec::from_iter(accessions.clone());

        Ok(Self {
            urs: parsed.urs().to_string(),
            taxid: parsed.taxid().try_into()?,
            short_urs: parsed.short(),
            deleted: String::from("N"),
            so_rna_type_tree: raw.so_tree().to_owned(),
            publication_count: raw.publication_count(),
            pre_summary: raw.precompute().into(),
            basic: base,
            qa_status: raw.qa_status().to_owned(),
            secondary: raw.r2dt().to_owned(),
            accessions: AccessionVec::from_iter(accessions.clone()),
            cross_references: accessions.into_iter().map(CrossReference::from).collect(),
            crs: raw.crs().to_owned().into_iter().collect(),
            overlaps: raw.feedback().to_owned().into_iter().collect(),
            go_annotations: raw.go_annotations().to_vec(),
            interacting_proteins: raw.interacting_proteins().to_vec(),
            interacting_rnas: raw.interacting_rnas().to_vec(),
            references,
            rfam_hits: raw.rfam_hits().to_owned().into_iter().collect(),
            orfs: raw.orfs().to_vec().into_iter().collect(),
            litsumm: raw.litsumm_summaries().to_vec(),
            editing_events: raw.editing_events().to_vec(),
        })
    }

    pub fn id(&self) -> &usize {
        &self.basic.id
    }

    pub fn urs_taxid(&self) -> &str {
        &self.basic.urs_taxid()
    }

    pub fn cross_references(&self) -> &Vec<CrossReference> {
        &self.cross_references
    }

    pub fn description(&self) -> &str {
        self.pre_summary.description()
    }
}
