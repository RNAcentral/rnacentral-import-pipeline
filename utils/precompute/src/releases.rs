use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::cmp::Ordering::{Less, Greater, Equal};

use serde::{Deserialize, Serialize};

use itertools::Itertools;
use sorted_iter::assume::*;
use sorted_iter::SortedPairIterator;

use anyhow::Result;

#[derive(Serialize, Deserialize, Debug)]
pub struct UrsEntry {
    urs: String,
    release: usize,
}

pub fn write_max(filename: &Path, output: &Path) -> Result<()> {
    let reader = rnc_utils::buf_reader(filename)?;
    let mut reader = csv::Reader::from_reader(reader);
    let out = rnc_utils::buf_writer(output)?;
    let mut writer = csv::Writer::from_writer(out);

    let results = reader
        .deserialize()
        .filter_map(Result::ok)
        .group_by(|e: &UrsEntry| e.urs.to_owned());

    for (_, entries) in &results {
        let max = entries.max_by(|l, r| l.release.cmp(&r.release));
        match max {
            Some(v) => writer.serialize(v)?,
            None => (),
        }
    }

    Ok(())
}

pub fn select_new(xrefs: &Path, known: &Path, output: &Path) -> Result<()> {
    let xref_buf = BufReader::new(File::open(xrefs)?);
    let mut xref_reader = csv::Reader::from_reader(xref_buf);
    let xref_records = xref_reader
        .deserialize()
        .filter_map(Result::ok)
        .map(|e: UrsEntry| (e.urs.to_owned(), e))
        .assume_sorted_by_key();

    let known_buf = BufReader::new(File::open(known)?);
    let mut known_reader = csv::Reader::from_reader(known_buf);
    let known_records = known_reader
        .deserialize()
        .filter_map(Result::ok)
        .map(|e: UrsEntry| (e.urs.to_owned(), e))
        .assume_sorted_by_key();

    let mut writer = csv::Writer::from_writer(File::create(output)?);
    let pairs = xref_records.outer_join(known_records);
    for (_key, (xref, pre)) in pairs {
        match (xref, pre) {
            (Some(x), Some(p)) => match x.release.cmp(&p.release) {
                Less => (),
                Equal => (),
                Greater => writer.write_record(&[x.urs])?,
            }
            (Some(x), None) => writer.write_record(&[x.urs])?,
            (None, Some(_)) => (),
            (None, None) => (),
        }
    }

    Ok(())
}
