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
    editing_events::EditingEvent,
    feedback::Feedback,
    go_annotation::GoAnnotation,
    go_flow_annotations::GoFlowLLMAnnotation,
    interacting_protein::InteractingProtein,
    interacting_rna::InteractingRna,
    litsumm::LitsummSummaries,
    orf::Orf,
    precompute::Precompute,
    publication_counts::PublicationCount,
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
    PublicationCount,
    LitsummSummaries,
    EditingEvents,
    GoFlowLLMAnnotations,
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
    publication_counts: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<PublicationCount>>,
    lit_summ: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<LitsummSummaries>>,
    go_flow_llm_annotations:
        StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<GoFlowLLMAnnotation>>,
    editing_events: StreamDeserializer<'de, IoRead<BufReader<File>>, Grouped<EditingEvent>>,
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
        let publication_counts = self.iterator_for(FileTypes::PublicationCount)?;
        let lit_summ = self.iterator_for(FileTypes::LitsummSummaries)?;
        let editing_events = self.iterator_for(FileTypes::EditingEvents)?;
        let go_flow_llm_annotations = self.iterator_for(FileTypes::GoFlowLLMAnnotations)?;
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
            publication_counts,
            lit_summ,
            editing_events,
            go_flow_llm_annotations,
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
            self.publication_counts.next(),
            self.lit_summ.next(),
            self.editing_events.next(),
            self.go_flow_llm_annotations.next(),
        );

        match current {
            (
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            ) => None,
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
                Some(Ok(Optional {
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
                Some(Ok(Optional {
                    id: id12,
                    data: publication_counts,
                })),
                Some(Ok(Multiple {
                    id: id13,
                    data: lit_summ,
                })),
                Some(Ok(Multiple {
                    id: id14,
                    data: editing_events,
                })),
                Some(Ok(Multiple {
                    id: id15,
                    data: goflow_llm_annotations,
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
                    || id1 != id12
                    || id1 != id13
                    || id1 != id14
                    || id1 != id15
                {
                    return Some(Err(Error::OutofSyncData(vec![
                        id1, id2, id3, id4, id5, id6, id7, id8, id9, id10, id11, id12, id13, id14,
                        id15,
                    ])));
                }

                let pre_so_type = precompute.so_rna_type();
                if !self.so_info.has_tree(&pre_so_type) {
                    return Some(Err(Error::MissingSoTrem(pre_so_type.to_string())));
                }
                let so_tree = self.so_info.tree(&pre_so_type).unwrap();

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
                    .publication_counts(publication_counts)
                    .litsumm_summaries(lit_summ)
                    .editing_events(editing_events)
                    .go_flow_llm_annotations(goflow_llm_annotations)
                    .so_tree(so_tree)
                    .build();

                Some(Ok(raw))
            },
            _ => Some(Err(Error::InvalidDataFormat)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    const SO_TREE: &str = r#"{"so_rna_type":"SO:0000655","so_term_tree":[["SO:0000655","ncRNA"]]}"#;

    fn base(id: usize) -> String {
        format!(
            r#"{{"Required":{{"id":{id},"data":{{"id":{id},"length":10,"md5":"abc","urs_taxid":"URS0000000001_9606"}}}}}}"#
        )
    }

    fn precompute(id: usize, so: &str) -> String {
        format!(
            r#"{{"Required":{{"id":{id},"data":{{"id":{id},"description":"d","rna_type":"ncRNA","has_coordinates":false,"so_rna_type":"{so}","databases":"ENA"}}}}}}"#
        )
    }

    fn multiple(id: usize) -> String {
        format!(r#"{{"Multiple":{{"id":{id},"data":[]}}}}"#)
    }

    fn optional(id: usize) -> String {
        format!(r#"{{"Optional":{{"id":{id},"data":null}}}}"#)
    }

    fn line_for(file: &FileTypes, base_id: usize, pre_id: usize, so: &str) -> String {
        match file {
            FileTypes::Base => base(base_id),
            FileTypes::Precompute => precompute(pre_id, so),
            FileTypes::QaStatus | FileTypes::R2dt | FileTypes::PublicationCount => optional(1),
            FileTypes::SoTermTree => SO_TREE.to_string(),
            _ => multiple(1),
        }
    }

    // Writes one file per type, named the way the workflow names them.
    fn write_files(dir: &Path, base_id: usize, pre_id: usize, so: &str) -> Vec<PathBuf> {
        FileTypes::iter()
            .map(|file| {
                let path = dir.join(format!("{}.json", file));
                fs::write(&path, line_for(&file, base_id, pre_id, so) + "\n").unwrap();
                path
            })
            .collect()
    }

    fn first(paths: Vec<PathBuf>) -> Result<Raw, Error> {
        FileJoinerBuilder::from_files(paths).unwrap().build().unwrap().next().unwrap()
    }

    #[test]
    fn file_names_map_to_their_type() {
        let cases = [
            ("base.json", FileTypes::Base),
            ("qa-status.json", FileTypes::QaStatus),
            ("rfam-hits.json", FileTypes::RfamHits),
            ("so-term-tree.json", FileTypes::SoTermTree),
            ("go-flow-llm-annotations.json", FileTypes::GoFlowLLMAnnotations),
        ];
        for (name, expected) in cases {
            assert_eq!(FileTypes::try_from(Path::new(name)).unwrap(), expected, "{}", name);
        }
    }

    #[test]
    fn an_unknown_file_name_is_rejected() {
        let err = FileTypes::try_from(Path::new("mystery.json")).unwrap_err();
        assert!(matches!(err, Error::UnknownFileType(_)));
    }

    #[test]
    fn every_file_type_must_be_given() {
        let dir = tempfile::tempdir().unwrap();
        let mut paths = write_files(dir.path(), 1, 1, "SO:0000655");
        paths.retain(|p| p.file_name().unwrap() != "orfs.json");
        let err = FileJoinerBuilder::from_files(paths).err().unwrap();
        assert!(matches!(err, Error::MustDefineMissingFile(FileTypes::Orfs)));
    }

    #[test]
    fn matching_ids_are_joined_into_one_raw_entry() {
        let dir = tempfile::tempdir().unwrap();
        let raw = first(write_files(dir.path(), 1, 1, "SO:0000655")).unwrap();
        assert_eq!(raw.id(), 1);
        assert_eq!(raw.urs_taxid(), "URS0000000001_9606");
        assert!(raw.crs().is_empty());
        assert!(raw.qa_status().is_none());
        assert_eq!(raw.so_tree().term_name().unwrap().to_string(), "ncRNA");
    }

    #[test]
    fn a_retired_so_term_is_looked_up_under_its_replacement() {
        let dir = tempfile::tempdir().unwrap();
        let raw = first(write_files(dir.path(), 1, 1, "SO:0001263")).unwrap();
        assert_eq!(raw.so_tree().term_name().unwrap().to_string(), "ncRNA");
    }

    #[test]
    fn files_that_disagree_on_the_id_are_reported_as_out_of_sync() {
        let dir = tempfile::tempdir().unwrap();
        let err = first(write_files(dir.path(), 1, 2, "SO:0000655")).unwrap_err();
        assert!(matches!(err, Error::OutofSyncData(ids) if ids[0] == 1 && ids[7] == 2));
    }

    #[test]
    fn a_missing_so_tree_is_reported_with_the_term() {
        let dir = tempfile::tempdir().unwrap();
        let err = first(write_files(dir.path(), 1, 1, "SO:0000374")).unwrap_err();
        assert!(matches!(err, Error::MissingSoTrem(term) if term == "SO:0000374"));
    }

    #[test]
    fn the_iterator_ends_when_every_file_ends() {
        let dir = tempfile::tempdir().unwrap();
        let paths = write_files(dir.path(), 1, 1, "SO:0000655");
        let mut joiner = FileJoinerBuilder::from_files(paths).unwrap().build().unwrap();
        assert!(joiner.next().unwrap().is_ok());
        assert!(matches!(joiner.next(), None));
    }

    #[test]
    fn a_file_with_extra_rows_is_an_invalid_format() {
        let dir = tempfile::tempdir().unwrap();
        let paths = write_files(dir.path(), 1, 1, "SO:0000655");
        let base_path = dir.path().join("base.json");
        fs::write(&base_path, base(1) + "\n" + &base(2) + "\n").unwrap();
        let mut joiner = FileJoinerBuilder::from_files(paths).unwrap().build().unwrap();
        assert!(joiner.next().unwrap().is_ok());
        assert!(matches!(joiner.next().unwrap(), Err(Error::InvalidDataFormat)));
    }
}
