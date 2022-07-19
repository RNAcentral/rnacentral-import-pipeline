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

use polars::prelude::*;

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

    let mut xref_records : DataFrame = CsvReader::from_path(xrefs)?.has_header(false).finish().unwrap();
    xref_records.rename("column_1", "id").ok();
    xref_records.rename("column_2", "upi").ok();
    xref_records.rename("column_3", "last").ok();
    let mut known_records : DataFrame = CsvReader::from_path(known)?.has_header(false).finish().unwrap();
    known_records.rename("column_1", "id").ok();
    known_records.rename("column_2", "upi").ok();
    known_records.rename("column_3", "last").ok();
    // Run groupby, sort and max on xref (because the DB doesn't have the memory to do it)
    xref_records = xref_records.groupby(["id", "upi"])?
                .select(["last"])
                .max()?
                .sort(["id"], false)
                .unwrap();

    // Join the frames on id, then filter to select those where xref > known (?)
    let mut selection = xref_records.join(&known_records, ["id", "upi"], ["id", "upi"], JoinType::Outer, None)?;
    let mask = selection.column("last_max")?.gt(selection.column("last")?)?;
    let mut selected_upis = selection.filter(&mask).unwrap()
                                 .select(["upi"])?
                                 .unique(None, UniqueKeepStrategy::First)?;


    let out_stream : File = File::create(output).unwrap();
    CsvWriter::new(out_stream)
        .has_header(false)
        .finish(&mut selected_upis);

    Ok(())
}
