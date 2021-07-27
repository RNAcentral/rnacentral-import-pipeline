use std::{
    collections::HashMap,
    fs::File,
    path::PathBuf,
    str::FromStr,
};

use serde::de::DeserializeOwned;
use strum_macros::{
    Display,
    EnumString,
};
use thiserror::Error;

use rnc_core::{
    grouper::Grouped::{
        self,
        Multiple,
        Optional,
        Required,
    },
    psql::JsonlIterator,
};

use super::{
    basic::Basic,
    crs::Crs,
    entry::Raw,
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
    so_tree::SoMapping,
};

#[derive(Debug, Error)]
pub enum Error {
    #[error("Missing required file {0}")]
    MustDefineMissingFile(FileTypes),

    #[error("Not all files are in required format")]
    InvalidDataFormat,

    #[error("Missing requested SO term {0}")]
    MissingSoTrem(String),

    #[error("Cannot find expected file, from {0:?}")]
    BadFileName(PathBuf),

    #[error("Could not match path {0:?} to a possible file type")]
    UnknownFileType(PathBuf),

    #[error("Data out of sync at: {0:?}")]
    OutofSyncData(Vec<usize>),

    #[error(transparent)]
    JsonError(#[from] rnc_core::psql::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

#[derive(Debug, Display, PartialEq, Eq, Hash, EnumString)]
pub enum FileTypes {
    #[strum(ascii_case_insensitive)]
    Base,

    #[strum(ascii_case_insensitive)]
    Crs,

    #[strum(ascii_case_insensitive)]
    Feedback,

    #[strum(ascii_case_insensitive)]
    GoAnnotations,

    #[strum(ascii_case_insensitive)]
    InteractingProteins,

    #[strum(ascii_case_insensitive)]
    InteractingRnas,

    #[strum(ascii_case_insensitive)]
    Orfs,

    #[strum(ascii_case_insensitive)]
    Precompute,

    #[strum(ascii_case_insensitive)]
    QaStatus,

    #[strum(ascii_case_insensitive)]
    R2dtHits,

    #[strum(ascii_case_insensitive)]
    RfamHits,

    #[strum(ascii_case_insensitive)]
    SoInfo,
}

pub struct FileJoiner {
    basic: JsonlIterator<File, Grouped<Basic>>,
    crs: JsonlIterator<File, Grouped<Crs>>,
    feedback: JsonlIterator<File, Grouped<Feedback>>,
    go_annotations: JsonlIterator<File, Grouped<GoAnnotation>>,
    interacting_proteins: JsonlIterator<File, Grouped<InteractingProtein>>,
    interacting_rnas: JsonlIterator<File, Grouped<InteractingRna>>,
    orfs: JsonlIterator<File, Grouped<Orf>>,
    precompute: JsonlIterator<File, Grouped<Precompute>>,
    qa_status: JsonlIterator<File, Grouped<QaStatus>>,
    r2dt_hits: JsonlIterator<File, Grouped<R2dt>>,
    rfam_hits: JsonlIterator<File, Grouped<RfamHit>>,
    so_info: SoMapping,
}

pub struct FileJoinerBuilder {
    paths: HashMap<FileTypes, PathBuf>,
}

impl Default for FileJoinerBuilder {
    fn default() -> Self {
        Self {
            paths: HashMap::default(),
        }
    }
}

impl FileJoinerBuilder {
    pub fn from_files(paths: Vec<PathBuf>) -> Result<Self, Error> {
        let mut builder = Self::default();
        for path in paths {
            let name = path
                .file_stem()
                .map(|s| s.to_str())
                .flatten()
                .ok_or_else(|| Error::BadFileName(path.clone()))?;
            let name = name.replace("-", "").replace("_", "");
            let file_type =
                FileTypes::from_str(&name).map_err(|_| Error::UnknownFileType(path.clone()))?;
            builder.file(file_type, path);
        }

        Ok(builder)
    }

    pub fn file(&mut self, file: FileTypes, path: PathBuf) -> &mut Self {
        self.paths.insert(file, path);
        self
    }

    fn path_for(&self, file: FileTypes) -> Result<&PathBuf, Error> {
        self.paths.get(&file).ok_or_else(|| Error::MustDefineMissingFile(file))
    }

    fn iterator_for<T: DeserializeOwned>(
        &self,
        file: FileTypes,
    ) -> Result<JsonlIterator<File, T>, Error> {
        let path = self.path_for(file)?;
        Ok(JsonlIterator::from_path(path)?)
    }

    pub fn build(&self) -> Result<FileJoiner, Error> {
        let basic = self.iterator_for(FileTypes::Base)?;
        let crs = self.iterator_for(FileTypes::Crs)?;
        let feedback = self.iterator_for(FileTypes::Feedback)?;
        let go_annotations = self.iterator_for(FileTypes::GoAnnotations)?;
        let interacting_proteins = self.iterator_for(FileTypes::InteractingProteins)?;
        let interacting_rnas = self.iterator_for(FileTypes::InteractingRnas)?;
        let orfs = self.iterator_for(FileTypes::Orfs)?;
        let precompute = self.iterator_for(FileTypes::Precompute)?;
        let qa_status = self.iterator_for(FileTypes::QaStatus)?;
        let r2dt_hits = self.iterator_for(FileTypes::R2dtHits)?;
        let rfam_hits = self.iterator_for(FileTypes::RfamHits)?;
        let so_info = so_tree::load(self.path_for(FileTypes::SoInfo)?)?;

        Ok(FileJoiner {
            basic,
            crs,
            feedback,
            go_annotations,
            interacting_proteins,
            interacting_rnas,
            orfs,
            precompute,
            qa_status,
            r2dt_hits,
            rfam_hits,
            so_info,
        })
    }
}

impl Iterator for FileJoiner {
    type Item = Result<Raw, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let current = (
            self.basic.next(),
            self.crs.next(),
            self.feedback.next(),
            self.go_annotations.next(),
            self.interacting_proteins.next(),
            self.interacting_rnas.next(),
            self.orfs.next(),
            self.precompute.next(),
            self.qa_status.next(),
            self.r2dt_hits.next(),
            self.rfam_hits.next(),
        );

        match current {
            (None, None, None, None, None, None, None, None, None, None, None) => None,
            (
                Some(Required {
                    id: id1,
                    data: base,
                }),
                Some(Multiple {
                    id: id2,
                    data: crs,
                }),
                Some(Multiple {
                    id: id3,
                    data: feedback,
                }),
                Some(Multiple {
                    id: id4,
                    data: go_annotations,
                }),
                Some(Multiple {
                    id: id5,
                    data: interacting_proteins,
                }),
                Some(Multiple {
                    id: id6,
                    data: interacting_rnas,
                }),
                Some(Multiple {
                    id: id7,
                    data: orfs,
                }),
                Some(Required {
                    id: id8,
                    data: precompute,
                }),
                Some(Required {
                    id: id9,
                    data: qa_status,
                }),
                Some(Optional {
                    id: id10,
                    data: r2dt,
                }),
                Some(Multiple {
                    id: id11,
                    data: rfam_hits,
                }),
            ) => {
                if id1 != id2
                    || id1 != id3
                    || id1 != id4
                    || id1 != id5
                    || id1 != id6
                    || id1 != id7
                    || id1 != id8
                    || id1 != id9
                    || id1 != id10
                    || id1 != id11
                {
                    return Some(Err(Error::OutofSyncData(vec![
                        id1, id2, id3, id4, id5, id6, id7, id8, id9, id10, id11,
                    ])));
                }

                let pre_so_type = precompute.so_rna_type();
                if !self.so_info.contains_key(pre_so_type) {
                    return Some(Err(Error::MissingSoTrem(pre_so_type.to_string())));
                }
                let so_tree = self.so_info[pre_so_type].clone();

                return Some(Ok(Raw {
                    id: id1,
                    base,
                    precompute,
                    qa_status,
                    crs,
                    feedback,
                    go_annotations,
                    interacting_proteins,
                    interacting_rnas,
                    r2dt,
                    rfam_hits,
                    orfs,
                    so_tree,
                }));
            },
            _ => Some(Err(Error::InvalidDataFormat)),
        }
    }
}