use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use serde::{Deserialize, Serialize};

use anyhow::Result;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoTreeEntry(String, String);

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SoTree(Vec<SoTreeEntry>);

pub type SoMapping = HashMap<String, SoTree>;

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
struct SoLine {
    so_rna_type: String,
    so_term_tree: SoTree,
}

impl SoTree {
    pub fn entries(&self) -> &[SoTreeEntry] {
        &self.0
    }

    /// Fetch the name of the SO term the tree is for.
    pub fn term_name(&self) -> Option<&str> {
        let terms = &self.0;
        if let Some(entry) = terms.last() {
            return Some(entry.name());
        }
        None
    }

    /// Fetch the SO id of the SO term the tree is for.
    pub fn term_id(&self) -> Option<&str> {
        let terms = &self.0;
        if let Some(entry) = terms.last() {
            return Some(entry.so_id());
        }
        None
    }

    /// Get a vector of all the names.
    pub fn names(&self) -> Vec<String> {
        let mut names = Vec::new();
        for entry in self.entries() {
            names.push(entry.name().to_string());
        }
        names
    }
}

impl SoTreeEntry {
    pub fn name(&self) -> &str {
        &self.1
    }

    pub fn so_id(&self) -> &str {
        &self.0
    }

    pub fn into_inner(self) -> (String, String) {
        (self.0, self.1)
    }
}

pub fn load(filename: &Path) -> Result<SoMapping> {
    let mut map: HashMap<String, SoTree> = HashMap::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line?;
        let entry: SoLine = serde_json::from_str(&line)?;
        map.insert(entry.so_rna_type, entry.so_term_tree);
    }

    Ok(map)
}
