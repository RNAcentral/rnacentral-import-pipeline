use std::{
    fs::File,
    io,
    io::{
        BufRead,
        BufReader,
    },
    path::Path,
};

use fnv::{
    FnvHashMap,
    FnvHashSet,
};

use thiserror::Error;

use crate::{
    urs::Urs,
    urs_taxid,
    urs_taxid::UrsTaxid,
};

#[derive(Error, Debug)]
pub enum Error {
    #[error("IO Error")]
    Io(#[from] io::Error),

    #[error("Incorrectly formatted URS_taxid")]
    UrsTaxidError(#[from] urs_taxid::Error),
}

pub struct UrsTaxidMapping {
    mapping: FnvHashMap<u64, FnvHashSet<u64>>,
}

impl Default for UrsTaxidMapping {
    fn default() -> Self {
        Self {
            mapping: FnvHashMap::default(),
        }
    }
}

impl UrsTaxidMapping {
    pub fn from_urs_file(path: &Path) -> Result<Self, Error> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let mut mapping: FnvHashMap<u64, FnvHashSet<u64>> = FnvHashMap::default();
        let mut buf = String::new();
        loop {
            match reader.read_line(&mut buf)? {
                0 => break,
                _ => {
                    let urs_taxid: UrsTaxid = buf.trim_end().parse()?;
                    let set = mapping.entry(urs_taxid.urs().into()).or_insert(FnvHashSet::default());
                    set.insert(urs_taxid.taxid());
                    buf.clear();
                },
            }
        }

        Ok(Self {
            mapping,
        })
    }

    pub fn urs_taxids(&self, urs: &Urs) -> Vec<UrsTaxid> {
        let id: u64 = urs.into();
        match self.mapping.get(&id) {
            None => Vec::with_capacity(0),
            Some(ts) => {
                let mut found = Vec::with_capacity(ts.len());
                for taxid in ts {
                    found.push(UrsTaxid::new(id, *taxid));
                }
                found
            },
        }
    }
}
