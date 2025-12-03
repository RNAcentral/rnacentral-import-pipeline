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

use crate::{
    fields::{
        SequenceEntry,
        SequenceFields,
        SoRnaTreeField,
    },
    search_xml::{
        SearchEntry,
        SearchValue,
    },
    sequences::{
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
    },
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
    qa_status: Option<QaStatus>,
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
            crs: raw.crs().iter().cloned().collect(),
            overlaps: raw.feedback().iter().cloned().collect(),
            go_annotations: raw.go_annotations().to_vec(),
            interacting_proteins: raw.interacting_proteins().to_vec(),
            interacting_rnas: raw.interacting_rnas().to_vec(),
            references,
            rfam_hits: raw.rfam_hits().iter().cloned().collect(),
            orfs: raw.orfs().iter().cloned().collect(),
            litsumm: raw.litsumm_summaries().to_vec(),
            editing_events: raw.editing_events().to_vec(),
        })
    }

    pub fn id(&self) -> &usize {
        &self.basic.id
    }

    pub fn urs_taxid(&self) -> &str {
        self.basic.urs_taxid()
    }

    pub fn cross_references(&self) -> &Vec<CrossReference> {
        &self.cross_references
    }

    pub fn description(&self) -> &str {
        self.pre_summary.description()
    }

    pub fn accessions(&self) -> &AccessionVec {
        &self.accessions
    }

    pub fn so_rna_type_tree(&self) -> &so_tree::SoTree {
        &self.so_rna_type_tree
    }

    pub fn has_litsumm(&self) -> bool {
        !self.litsumm.is_empty()
    }

    pub fn has_lit_scan(&self) -> bool {
        self.publication_count > 0
    }

    pub fn has_go_annotations(&self) -> bool {
        !self.go_annotations.is_empty()
    }

    pub fn has_secondary_structure(&self) -> bool {
        self.secondary.is_some()
    }

    pub fn has_qc_warning(&self) -> bool {
        if let Some(qa) = &self.qa_status {
            return qa.has_issue();
        }
        false
    }

    pub fn precompute_summary(&self) -> &PrecomputeSummary {
        &self.pre_summary
    }
}

impl SearchEntry<SequenceEntry> for Normalized {
    fn id(&self) -> &str {
        self.urs_taxid()
    }

    fn name(&self) -> &str {
        self.urs_taxid()
    }

    fn description(&self) -> &str {
        self.description()
    }

    fn taxid(&self) -> usize {
        self.taxid
    }

    fn field_values(&self, field: &SequenceFields) -> SearchValue {
        match field {
            SequenceFields::ShortUrs => self.short_urs.clone().into(),
            SequenceFields::Length => self.basic.length.into(),
            SequenceFields::Organelle => todo!(),
            SequenceFields::ExpertDb => self.pre_summary.databases().into(),
            SequenceFields::CommonName => todo!(),
            SequenceFields::Function => todo!(),
            SequenceFields::Gene => todo!(),
            SequenceFields::GeneSynonym => todo!(),
            SequenceFields::InsdcRNAType => todo!(),
            SequenceFields::Product => todo!(),
            SequenceFields::HasGenomicCoordinates => self.pre_summary.has_coordinates().into(),
            SequenceFields::Md5 => self.basic.md5.clone().into(),
            SequenceFields::Author => todo!(),
            SequenceFields::Journal => todo!(),
            SequenceFields::InsdcSubmission => todo!(),
            SequenceFields::PubTitle => todo!(),
            SequenceFields::PubId => todo!(),
            SequenceFields::PopularSpecies => todo!(),
            SequenceFields::Boost => todo!(),
            SequenceFields::LocusTag => todo!(),
            SequenceFields::StandardName => todo!(),
            SequenceFields::RfamFamilyName => todo!(),
            SequenceFields::RfamId => todo!(),
            SequenceFields::RfamClan => todo!(),
            SequenceFields::QcWarning => todo!(),
            SequenceFields::QcWarningFound => todo!(),
            SequenceFields::TaxString => todo!(),
            SequenceFields::InvovledIn => todo!(),
            SequenceFields::PartOf => todo!(),
            SequenceFields::Enables => todo!(),
            SequenceFields::ContributesTo => todo!(),
            SequenceFields::ColocalizesWith => todo!(),
            SequenceFields::HasGoAnnotations => todo!(),
            SequenceFields::GoAnnotationSource => todo!(),
            SequenceFields::HasInteractingProteins => todo!(),
            SequenceFields::InteractingProtein => todo!(),
            SequenceFields::HasInteractingRnas => todo!(),
            SequenceFields::HasConservedStructure => todo!(),
            SequenceFields::ConservedStructure => todo!(),
            SequenceFields::OverlapsWith => todo!(),
            SequenceFields::NoOverlapsWith => todo!(),
            SequenceFields::HasSecondaryStructure => todo!(),
            SequenceFields::SecondaryStructureModel => todo!(),
            SequenceFields::SecondaryStructureSource => todo!(),
            SequenceFields::PdbidEntityid => todo!(),
            SequenceFields::Disease => todo!(),
            SequenceFields::Url => todo!(),
            SequenceFields::OrfSource => todo!(),
            SequenceFields::HasLitScan => todo!(),
            SequenceFields::HasLitsumm => todo!(),
            SequenceFields::HasEditingEvent => todo!(),
            SequenceFields::EditChromosome => todo!(),
            SequenceFields::EditLocations => todo!(),
            SequenceFields::EditRepeatType => todo!(),
            SequenceFields::SoRnaTypeName => todo!(),
            SequenceFields::SoRnaType => todo!(),
        }
    }

    fn tree_field_values(&self, field: &SoRnaTreeField) -> (String, Vec<String>) {
        match field {
            SoRnaTreeField::SoRnaType => todo!(),
        }
    }

    fn cross_references(&self) -> impl IntoIterator<Item = CrossReference> {
        self.cross_references.clone()
    }
}
