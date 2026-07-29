use core::{
    convert::{
        Into,
        TryFrom,
    },
    num::{
        NonZeroU16,
        NonZeroU32,
    },
};
use std::{
    collections::HashMap,
    fmt,
    fs::File,
    io::{
        BufRead,
        BufReader,
    },
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;
use thiserror::Error;

use crate::search_xml::AsSearchValue;

#[derive(Error, Debug)]
pub enum SoIdError {
    #[error("SO ID must start with 'SO:', got: {0}")]
    MissingPrefix(String),

    #[error("SO ID must have exactly 7 digits after 'SO:', got {0} digits: `{1}")]
    InvalidLength(usize, String),

    #[error("Failed to parse SO ID number: {0}")]
    ParseError(#[from] std::num::ParseIntError),

    #[error("SO ID number cannot be zero")]
    ZeroValue,
}

/// This is a wrapper around SO ids, like `SO:0000374`. This is used to track and esnure
/// that the ids are used where they are meant to and names where they should be.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Hash)]
#[serde(try_from = "String", into = "String")]
pub struct SoId(NonZeroU32);

/// This wraps SO names.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoName(String);

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoTreeEntry(SoId, SoName);

/// An `SoTree` represents the
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoTree(Vec<SoTreeEntry>);

pub struct SoMapping(HashMap<SoId, SoTree>);

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
struct SoLine {
    so_rna_type: SoId,
    so_term_tree: SoTree,
}

impl fmt::Display for SoId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SO:{:07}", self.0.get())
    }
}

impl TryFrom<String> for SoId {
    type Error = SoIdError;

    fn try_from(value: String) -> core::result::Result<Self, Self::Error> {
        // Check if string starts with "SO:"
        let suffix =
            value.strip_prefix("SO:").ok_or_else(|| SoIdError::MissingPrefix(value.clone()))?;

        // Check if suffix is exactly 7 digits
        if suffix.len() != 7 {
            return Err(SoIdError::InvalidLength(suffix.len(), suffix.to_string()));
        }

        // Parse the numeric part (ParseIntError automatically converted via #[from])
        let num = suffix.parse::<u32>()?;

        // Create NonZeroU16
        let non_zero = NonZeroU32::new(num).ok_or(SoIdError::ZeroValue)?;

        Ok(SoId(non_zero))
    }
}

impl TryFrom<&str> for SoId {
    type Error = SoIdError;

    fn try_from(value: &str) -> std::result::Result<Self, Self::Error> {
        Self::try_from(value.to_string())
    }
}

impl Into<String> for SoId {
    fn into(self) -> String {
        format!("SO:{:07}", self.0.get())
    }
}

impl AsRef<str> for SoName {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Into<String> for SoName {
    fn into(self) -> String {
        self.0
    }
}

impl fmt::Display for SoName {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl AsSearchValue for SoId {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl AsSearchValue for &SoId {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl AsSearchValue for SoName {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl AsSearchValue for &SoName {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl SoMapping {
    pub fn has_tree(&self, so_id: &SoId) -> bool {
        self.0.contains_key(so_id)
    }

    pub fn tree(&self, so_id: &SoId) -> Option<SoTree> {
        self.0.get(so_id).cloned()
    }
}

impl SoTree {
    pub fn entries(&self) -> &[SoTreeEntry] {
        &self.0
    }

    /// Fetch the name of the SO term the tree is for.
    pub fn term_name(&self) -> Option<&SoName> {
        let terms = &self.0;
        if let Some(entry) = terms.last() {
            return Some(entry.name());
        }
        None
    }

    /// Fetch the SO id of the SO term the tree is for.
    pub fn term_id(&self) -> Option<&SoId> {
        let terms = &self.0;
        if let Some(entry) = terms.last() {
            return Some(entry.so_id());
        }
        None
    }

    /// Get a vector of all the names.
    pub fn names(&self) -> Vec<&SoName> {
        let mut names = Vec::new();
        for entry in self.entries() {
            names.push(entry.name());
        }
        names
    }
}

impl SoTreeEntry {
    pub fn name(&self) -> &SoName {
        &self.1
    }

    pub fn so_id(&self) -> &SoId {
        &self.0
    }

    pub fn into_inner(self) -> (SoId, SoName) {
        (self.0, self.1)
    }
}

/// Load a flie of SO information. The file should be a JSONL (one JSON object per line)
/// formatted file, where each line is an SoLine.
pub fn load(filename: &Path) -> Result<SoMapping> {
    let mut map: HashMap<SoId, SoTree> = HashMap::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let entry: SoLine = serde_json::from_str(&line)?;
        map.insert(entry.so_rna_type, entry.so_term_tree);
    }

    Ok(SoMapping(map))
}
