use std::{
    collections::HashMap,
    convert::TryFrom,
    fs::File,
    io::BufReader,
    path::{
        Path,
        PathBuf,
    },
    str::FromStr,
};

use serde::de::DeserializeOwned;
use serde_json::{
    de::IoRead,
    Deserializer,
    StreamDeserializer,
};
use strum::IntoEnumIterator;
use strum_macros::{
    Display,
    EnumIter,
    EnumString,
};
use thiserror::Error;

use rnc_core::grouper::Grouped::{
    self,
    Multiple,
    Optional,
    Required,
};

use super::{
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
    raw::Raw,
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
    IoError(#[from] std::io::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

#[derive(Debug, Display, PartialEq, Eq, Hash, EnumString, EnumIter)]
#[strum(ascii_case_insensitive, serialize_all = "kebab-case")]
pub enum FileTypes {
    Base,
    Crs,
    Feedback,
    GoAnnotations,
    InteractingProteins,
    InteractingRnas,
    Orfs,
    Precompute,
    QaStatus,
    R2dt,
    RfamHits,
    SoTermTree,
}

pub struct FileJoiner<'de> {
    basic: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<Basic>>,
    crs: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<Crs>>,
    feedback: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<Feedback>>,
    go_annotations: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<GoAnnotation>>,
    interacting_proteins:
        StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<InteractingProtein>>,
    interacting_rnas: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<InteractingRna>>,
    orfs: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<Orf>>,
    precompute: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<Precompute>>,
    qa_status: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<QaStatus>>,
    r2dt_hits: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<R2dt>>,
    rfam_hits: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<RfamHit>>,
    so_info: SoMapping,
}

pub struct FileJoinerBuilder {
    paths: HashMap<FileTypes, PathBuf>,
}

impl TryFrom<&Path> for FileTypes {
    type Error = Error;

    fn try_from(path: &Path) -> Result<Self, Self::Error> {
        let name = path
            .file_stem()
            .map(|s| s.to_str())
            .flatten()
            .ok_or_else(|| Error::BadFileName(PathBuf::from(path)))?;
        Self::from_str(&name).map_err(|_| Error::UnknownFileType(PathBuf::from(path)))
    }
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
            let file_type = FileTypes::try_from(path.as_ref())?;
            builder.file(file_type, path);
        }

        for file_type in FileTypes::iter() {
            if !builder.is_set(&file_type) {
                return Err(Error::MustDefineMissingFile(file_type));
            }
        }

        Ok(builder)
    }

    pub fn is_set(&mut self, file: &FileTypes) -> bool {
        self.paths.contains_key(file)
    }

    pub fn file(&mut self, file: FileTypes, path: PathBuf) -> &mut Self {
        self.paths.insert(file, path);
        self
    }

    fn path_for(&self, file: FileTypes) -> Result<&PathBuf, Error> {
        self.paths.get(&file).ok_or_else(|| Error::MustDefineMissingFile(file))
    }

    fn iterator_for<'de, T>(
        &self,
        file: FileTypes,
    ) -> Result<StreamDeserializer<'de, IoRead<BufReader<File>>, T>, Error>
    where
        T: DeserializeOwned,
    {
        let path = self.path_for(file)?;
        let reader = BufReader::new(File::open(&path)?);
        let stream = Deserializer::from_reader(reader).into_iter::<T>();
        Ok(stream)
    }

    pub fn build<'de>(&self) -> Result<FileJoiner<'de>, Error> {
        let basic = self.iterator_for(FileTypes::Base)?;
        let crs = self.iterator_for(FileTypes::Crs)?;
        let feedback = self.iterator_for(FileTypes::Feedback)?;
        let go_annotations = self.iterator_for(FileTypes::GoAnnotations)?;
        let interacting_proteins = self.iterator_for(FileTypes::InteractingProteins)?;
        let interacting_rnas = self.iterator_for(FileTypes::InteractingRnas)?;
        let orfs = self.iterator_for(FileTypes::Orfs)?;
        let precompute = self.iterator_for(FileTypes::Precompute)?;
        let qa_status = self.iterator_for(FileTypes::QaStatus)?;
        let r2dt_hits = self.iterator_for(FileTypes::R2dt)?;
        let rfam_hits = self.iterator_for(FileTypes::RfamHits)?;
        let so_info = so_tree::load(self.path_for(FileTypes::SoTermTree)?)?;

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

impl<'de> Iterator for FileJoiner<'de> {
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
                Some(Ok(Required {
                    id: id1,
                    data: base,
                })),
                Some(Ok(Multiple {
                    id: id2,
                    data: crs,
                })),
                Some(Ok(Multiple {
                    id: id3,
                    data: feedback,
                })),
                Some(Ok(Multiple {
                    id: id4,
                    data: go_annotations,
                })),
                Some(Ok(Multiple {
                    id: id5,
                    data: interacting_proteins,
                })),
                Some(Ok(Multiple {
                    id: id6,
                    data: interacting_rnas,
                })),
                Some(Ok(Multiple {
                    id: id7,
                    data: orfs,
                })),
                Some(Ok(Required {
                    id: id8,
                    data: precompute,
                })),
                Some(Ok(Required {
                    id: id9,
                    data: qa_status,
                })),
                Some(Ok(Optional {
                    id: id10,
                    data: r2dt,
                })),
                Some(Ok(Multiple {
                    id: id11,
                    data: rfam_hits,
                })),
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

                let raw = Raw::builder()
                    .id(id1)
                    .base(base)
                    .precompute(precompute)
                    .qa_status(qa_status)
                    .crs(crs)
                    .feedback(feedback)
                    .go_annotations(go_annotations)
                    .interacting_proteins(interacting_proteins)
                    .interacting_rnas(interacting_rnas)
                    .r2dt(r2dt)
                    .rfam_hits(rfam_hits)
                    .orfs(orfs)
                    .so_tree(so_tree)
                    .build();

                Some(Ok(raw))
            },
            _ => Some(Err(Error::InvalidDataFormat)),
        }
    }
}
