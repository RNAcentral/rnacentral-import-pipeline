use std::iter::FromIterator;
use std::convert::TryInto;

use serde::{
    Deserialize,
    Serialize,
};

use thiserror::Error;

use rnc_core::{
    urs,
    urs_taxid,
};

use crate::normalize::{
    accession::{
        AccessionVec,
        CrossReference,
        RawAccession,
        ReferenceVec,
    },
    basic::Basic,
    crs::{
        Crs,
        CrsVec,
    },
    feedback::{
        Feedback,
        FeedbackVec,
    },
    go_annotation::GoAnnotation,
    interacting_protein::InteractingProtein,
    interacting_rna::InteractingRna,
    precompute::{
        Precompute,
        PrecomputeSummary,
    },
    qa_status::QaStatus,
    r2dt::R2dt,
    rfam_hit::{
        RfamHit,
        RfamHitVec,
    },
    so_tree,
};

#[derive(Error, Debug)]
pub enum NormalizationError {
    #[error("Could not parse {0}")]
    UrsParsingError(#[from] urs::Error),

    #[error("Could not parse {0}")]
    UrsTaxidParsingError(#[from] urs_taxid::Error),
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Raw {
    pub id: usize,
    pub base: Basic,
    pub crs: Vec<Crs>,
    pub feedback: Vec<Feedback>,
    pub go_annotations: Vec<GoAnnotation>,
    pub interacting_proteins: Vec<InteractingProtein>,
    pub interacting_rnas: Vec<InteractingRna>,
    pub precompute: Precompute,
    pub qa_status: QaStatus,
    pub r2dt: Option<R2dt>,
    pub rfam_hits: Vec<RfamHit>,
    pub so_tree: so_tree::SoTree,
}

#[derive(Debug, PartialEq, Serialize, Deserialize)]
pub struct Normalized {
    urs: String,
    taxid: usize,
    urs_taxid: String,
    short_urs: String,
    deleted: String,

    so_rna_type_tree: so_tree::SoTree,

    #[serde(flatten)]
    pre_summary: PrecomputeSummary,

    #[serde(flatten)]
    basic: Basic,

    // #[serde(flatten)]
    // dates: Dates,
    qa_status: QaStatus,

    secondary: Option<R2dt>,

    #[serde(flatten)]
    accessions: AccessionVec,
    cross_references: Vec<CrossReference>,
    crs: CrsVec,
    overlaps: FeedbackVec,
    go_annotations: Vec<GoAnnotation>,

    interacting_proteins: Vec<InteractingProtein>,

    interacting_rnas: Vec<InteractingRna>,

    #[serde(flatten)]
    references: ReferenceVec,

    #[serde(flatten)]
    rfam_hits: RfamHitVec,
}

impl Raw {
    pub fn urs_taxid(&self) -> String {
        return self.base.urs_taxid.to_owned();
    }

    pub fn urs(&self) -> Result<String, urs_taxid::Error> {
        let ut: urs_taxid::UrsTaxid = self.urs_taxid().parse()?;
        let urs: urs::Urs = ut.into();
        Ok(urs.to_string())
    }

    pub fn taxid(&self) -> Result<u64, urs_taxid::Error> {
        let ut: urs_taxid::UrsTaxid = self.urs_taxid().parse()?;
        Ok(ut.taxid())
    }

    pub fn short_urs(&self) -> Result<String, urs_taxid::Error> {
        let ut: urs_taxid::UrsTaxid = self.urs_taxid().parse()?;
        let urs: urs::Urs = ut.into();
        Ok(urs.short_urs())
    }

    pub fn short_urs_taxid(&self) -> Result<String, urs_taxid::Error> {
        let ut: urs_taxid::UrsTaxid = self.urs_taxid().parse()?;
        Ok(ut.short())
    }
}

impl Normalized {
    pub fn new(raw: Raw, accessions: Vec<RawAccession>) -> Result<Self, NormalizationError> {
        let pre_summary = PrecomputeSummary::from(raw.precompute);
        let base = raw.base.clone();

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
            urs_taxid,
            urs: parsed.urs().to_string(),
            taxid: parsed.taxid().try_into().unwrap(),
            short_urs: parsed.short(),
            deleted: String::from("N"),
            so_rna_type_tree: raw.so_tree,
            pre_summary,
            basic: base,
            qa_status: raw.qa_status,
            secondary: raw.r2dt,
            accessions: AccessionVec::from_iter(accessions.clone()),
            cross_references: accessions.into_iter().map(CrossReference::from).collect(),
            crs: CrsVec::from_iter(raw.crs.clone()),
            overlaps: FeedbackVec::from_iter(raw.feedback.clone()),
            go_annotations: raw.go_annotations.clone(),
            interacting_proteins: raw.interacting_proteins,
            interacting_rnas: raw.interacting_rnas,
            references,
            rfam_hits: RfamHitVec::from_iter(raw.rfam_hits.clone()),
        })
    }
}
