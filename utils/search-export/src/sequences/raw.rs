use serde::{
    Deserialize,
    Serialize,
};

use typed_builder::TypedBuilder;

use rnc_core::{
    urs,
    urs_taxid,
};

use crate::sequences::{
    basic::Basic,
    crs::Crs,
    feedback::Feedback,
    go_annotation::GoAnnotation,
    interacting_protein::InteractingProtein,
    interacting_rna::InteractingRna,
    orf::Orf,
    precompute::Precompute,
    qa_status::QaStatus,
    r2dt::R2dt,
    rfam_hit::RfamHit,
    so_tree,
};

#[derive(TypedBuilder, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Raw {
    id: usize,
    base: Basic,
    crs: Vec<Crs>,
    feedback: Vec<Feedback>,
    go_annotations: Vec<GoAnnotation>,
    interacting_proteins: Vec<InteractingProtein>,
    interacting_rnas: Vec<InteractingRna>,
    precompute: Precompute,
    qa_status: QaStatus,
    r2dt: Option<R2dt>,
    rfam_hits: Vec<RfamHit>,
    orfs: Vec<Orf>,
    so_tree: so_tree::SoTree,
}

impl Raw {
    pub fn id(&self) -> usize {
        self.id
    }

    pub fn urs_taxid(&self) -> &String {
        &self.base.urs_taxid
    }

    pub fn crs(&self) -> &[Crs] {
        &self.crs
    }

    pub fn feedback(&self) -> &[Feedback] {
        &self.feedback
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

    /// Get a reference to the raw's go annotations.
    pub fn go_annotations(&self) -> &[GoAnnotation] {
        self.go_annotations.as_slice()
    }

    /// Get a reference to the raw's interacting proteins.
    pub fn interacting_proteins(&self) -> &[InteractingProtein] {
        self.interacting_proteins.as_slice()
    }

    /// Get a reference to the raw's interacting rnas.
    pub fn interacting_rnas(&self) -> &[InteractingRna] {
        self.interacting_rnas.as_slice()
    }

    /// Get a reference to the raw's precompute.
    pub fn precompute(&self) -> &Precompute {
        &self.precompute
    }

    /// Get a reference to the raw's qa status.
    pub fn qa_status(&self) -> &QaStatus {
        &self.qa_status
    }

    /// Get a reference to the raw's r2dt.
    pub fn r2dt(&self) -> &Option<R2dt> {
        &self.r2dt
    }

    /// Get a reference to the raw's rfam hits.
    pub fn rfam_hits(&self) -> &[RfamHit] {
        self.rfam_hits.as_slice()
    }

    /// Get a reference to the raw's orfs.
    pub fn orfs(&self) -> &[Orf] {
        self.orfs.as_slice()
    }

    /// Get a reference to the raw's so tree.
    pub fn so_tree(&self) -> &so_tree::SoTree {
        &self.so_tree
    }

    /// Get a reference to the raw's base.
    pub fn base(&self) -> &Basic {
        &self.base
    }
}
