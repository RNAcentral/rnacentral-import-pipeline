use std::{
    cmp::Ordering::{
        Equal,
        Greater,
        Less,
    },
    fs::File,
    io::BufReader,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use itertools::Itertools;
use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use anyhow::{
    anyhow,
    Result,
};

#[derive(Serialize, Deserialize, Debug)]
pub struct UrsEntry {
    id: usize,
    urs: String,
    release: usize,
}

fn entries(path: &Path) -> Result<impl Iterator<Item = UrsEntry>> {
    let file = File::open(path)?;
    let buf = BufReader::new(file);
    let reader = csv::ReaderBuilder::new().has_headers(false).from_reader(buf);

    let records = reader.into_deserialize().filter_map(Result::ok);

    Ok(records)
}

pub fn write_max(filename: &Path, output: &Path) -> Result<()> {
    let out = rnc_utils::buf_writer(output)?;
    let mut writer = csv::Writer::from_writer(out);

    let results = entries(filename)?.group_by(|e: &UrsEntry| e.id);

    for (_, entries) in &results {
        let max = entries.max_by(|l, r| l.release.cmp(&r.release));
        match max {
            Some(v) => writer.serialize(v)?,
            None => (),
        }
    }
    writer.flush()?;

    Ok(())
}

pub fn select_new(xrefs: &Path, known: &Path, output: &Path) -> Result<()> {
    let xref_records = entries(xrefs)?.map(|e: UrsEntry| (e.id, e)).assume_sorted_by_key();
    let known_records = entries(known)?.map(|e: UrsEntry| (e.id, e)).assume_sorted_by_key();

    let mut writer = csv::Writer::from_writer(File::create(output)?);
    let pairs = xref_records.outer_join(known_records);
    for (_key, (xref, pre)) in pairs {
        match (xref, pre) {
            (Some(x), Some(p)) => match x.release.cmp(&p.release) {
                Less => Err(anyhow!(
                    "This should never happen, too small release for {:?} vs {:?}",
                    &x,
                    &p
                ))?,
                Equal => (),
                Greater => writer.write_record(&[x.urs])?,
            },
            (Some(x), None) => writer.write_record(&[x.urs])?,
            (None, Some(_)) => (),
            (None, None) => (),
        }
    }
    writer.flush()?;

    Ok(())
}
