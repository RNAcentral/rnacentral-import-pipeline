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
use serde_json::error;
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

/// Selects new UPIs from a set of cross-reference records and a set of known records.
///
/// # Description
///
/// Used to select wqhich UPIs are forwarded to the precompute workflow, based on when
/// they were last updated. Xref is updated as a sequence is imported, so the value of
/// last in there is the last time we updated the UPI. known comes from the existing
/// precompute table. If the last value in known is the same as the last value in xref,
/// then we don't need to update the UPI. If the last value in known is less than the last
/// value in xref, then we need to update the UPI. If the last value in known is greater
/// than the last value in xref, then something has gone wrong.
///
/// This function reads two CSV files containing cross-reference records and known
/// records, respectively. It then joins the two data frames on the `id` and `upi`
/// columns, filters the resulting data frame to only include rows where the `last_xref`
/// column is greater than the `last_precompute` column, selects only the `upi` column,
/// and removes any duplicate rows based on the `upi` column. The resulting data frame is
/// then written to a new CSV file specified by the `output` argument.
///
/// # Arguments
///
/// * `xrefs` - A `Path` to the CSV file containing cross-reference records.
/// * `known` - A `Path` to the CSV file containing known records.
/// * `output` - A `Path` to the CSV file to write the selected UPIs to.
///
/// # Errors
///
/// This function returns an error if there is an issue reading or writing the CSV files,
/// or if there is an issue with the data frames or data frame operations.
///
/// # Examples
///
/// ```
/// use std::path::Path;
/// use anyhow::Result;
/// use polars::prelude::*;
///
/// fn main() -> Result<()> {
///     let xrefs = Path::new("xrefs.csv");
///     let known = Path::new("known.csv");
///     let output = Path::new("selected_upis.csv");
///     select_new(&xrefs, &known, &output)?;
///     Ok(())
/// }
/// ```
pub fn select_new(xrefs: &Path, known: &Path, output: &Path, streaming: bool) -> Result<()> {
    // println!("{:?} {:?} {:?}", xrefs, known, streaming);
    // let known_path = known.to_str().unwrap().to_owned();
    // let xrefs_path = xrefs.to_str().unwrap().to_owned();
    // let xref_records: LazyFrame = LazyCsvReader::new(xrefs_path)
    //     .has_header(false)
    //     .low_memory(streaming)
    //     .finish()?
    //     .rename(vec!["column_1", "column_2", "column_3"], vec!["id", "upi", "last"])
    //     .group_by(["upi"])
    //     .agg([col("last").max().alias("last"), col("id").first().alias("id")])
    //     .sort("id", Default::default());


    // let known_records: LazyFrame = LazyCsvReader::new(known_path)
    //     .has_header(false)
    //     .low_memory(streaming)
    //     .finish()
    //     .unwrap()
    //     .rename(vec!["column_1", "column_2", "column_3"], vec!["id", "upi", "last"])
    //     .group_by(["upi"])
    //     .agg([col("last").max().alias("last"), col("id").first().alias("id")])
    //     .sort("id", Default::default());


    // let selection: LazyFrame = xref_records
    //     .join(
    //         known_records,
    //         [col("upi")],
    //         [col("upi")],
    //         JoinArgs::new(JoinType::Outer),
    //     )
    //     .rename(vec!["last", "last_right"], vec!["last_xref", "last_precompute"])
    //     .with_columns([
    //         col("last_xref").gt(col("last_precompute")).alias("selected"),
    //         col("last_precompute").gt(col("last_xref")).alias("error_state"),
    //     ]);
    //     // .select([col("upi"), col("selected"), col("error_state")]);

    // let check: LazyFrame = selection.clone();

    // // // check we are not in a catastrophic error state - precompute should never be newer than
    // // // xref
    // let selected_urs = selection.filter(col("selected").eq(true)).with_streaming(streaming).collect()?;
    // let error_urs = check.filter(col("error_state").eq(true)).with_streaming(streaming).collect()?;
    // if error_urs.height() > 0 {
    //     return Err(anyhow!("Precompute newer than xref for these UPIs: {:?}", error_urs));
    // }
    // println!("{:?}", selected_urs);

    // let mut selected_upis =
    //     selected_urs.select(["upi"])?.unique(None, UniqueKeepStrategy::First, None)?;

    // let out_stream: File = File::create(output).unwrap();
    // CsvWriter::new(out_stream).has_header(false).finish(&mut selected_upis)?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{
        distributions::Alphanumeric,
        Rng,
    };
    use std::io::Cursor; // 0.8

    fn get_random_fname() -> String {
        let rand_string: String =
            rand::thread_rng().sample_iter(&Alphanumeric).take(30).map(char::from).collect();
        let fname = format!("{}.csv", rand_string);
        fname
    }

    #[test]
    fn test_select_new() -> Result<()> {
        // Create test data frames
        let xref_data = "1,upi1,601\n\
                         2,upi2,600\n\
                         3,upi3,602\n\
                         4,upi4,600\n\
                         5,upi5,603\n\
                         1,upi1,604\n\
                         2,upi2,600\n\
                         3,upi3,605\n\
                         4,upi4,600\n\
                         5,upi5,606\n";
        let known_data = "1,upi1,600\n\
                          2,upi2,600\n\
                          3,upi3,600\n\
                          4,upi4,600\n\
                          5,upi5,600\n\
                          6,upi1,500\n\
                          7,upi2,600\n\
                          8,upi3,500\n\
                          9,upi4,600\n\
                          10,upi5,500\n";
        let xref_reader = CsvReader::new(Cursor::new(xref_data)).has_header(false);
        let known_reader = CsvReader::new(Cursor::new(known_data)).has_header(false);
        let mut xref_records = xref_reader.finish()?;
        let mut known_records = known_reader.finish()?;

        let xref_fname = get_random_fname();
        let known_fname = get_random_fname();

        // write the test data frames to files
        let mut xref_path = File::create(&xref_fname)?;
        let mut known_path = File::create(&known_fname)?;
        let xref_writer = CsvWriter::new(&mut xref_path);
        let known_writer = CsvWriter::new(&mut known_path);
        xref_writer.has_header(false).finish(&mut xref_records)?;
        known_writer.has_header(false).finish(&mut known_records)?;

        // Call the function being tested
        let output_fname = get_random_fname();
        let output = Path::new(&output_fname);
        select_new(Path::new(&xref_fname), Path::new(&known_fname), &output, false)?;

        // Read the output file and check its contents
        let output_reader = CsvReader::from_path(&output)?.has_header(false);
        let output_records = output_reader.finish()?;
        let expected_data = "upi1\n\
                             upi3\n\
                             upi5\n";
        let expected_reader = CsvReader::new(Cursor::new(expected_data)).has_header(false);
        let expected_records = expected_reader.finish()?;
        assert_eq!(output_records.sort(["column_1"], false, false).unwrap(), expected_records);

        // Clean up the output file
        std::fs::remove_file(&output)?;
        std::fs::remove_file(Path::new(&xref_fname))?;
        std::fs::remove_file(Path::new(&known_fname))?;

        Ok(())
    }

    #[test]
    fn test_select_new_known_newer() -> Result<()> {
        // Create test data frames
        // Create test data frames
        let xref_data = "1,upi1,600\n\
                         2,upi2,600\n\
                         3,upi3,600\n\
                         4,upi4,601\n\
                         5,upi5,600\n\
                         1,upi1,600\n\
                         2,upi2,600\n\
                         3,upi3,600\n\
                         4,upi4,501\n\
                         5,upi5,600\n";
        let known_data = "1,upi1,600\n\
                          2,upi2,600\n\
                          3,upi3,600\n\
                          4,upi4,602\n\
                          5,upi5,600\n\
                          1,upi1,600\n\
                          2,upi2,600\n\
                          3,upi3,600\n\
                          4,upi4,703\n\
                          5,upi5,600\n";
        let xref_reader = CsvReader::new(Cursor::new(xref_data)).has_header(false);
        let known_reader = CsvReader::new(Cursor::new(known_data)).has_header(false);
        let mut xref_records = xref_reader.finish()?;
        let mut known_records = known_reader.finish()?;

        let xref_fname = get_random_fname();
        let known_fname = get_random_fname();
        // write the test data frames to files
        let mut xref_path = File::create(&xref_fname)?;
        let mut known_path = File::create(&known_fname)?;
        let xref_writer = CsvWriter::new(&mut xref_path);
        let known_writer = CsvWriter::new(&mut known_path);
        xref_writer.has_header(false).finish(&mut xref_records)?;
        known_writer.has_header(false).finish(&mut known_records)?;

        // Call the function being tested
        let output_fname = get_random_fname();
        let output = Path::new(&output_fname);
        let result = select_new(Path::new(&xref_fname), Path::new(&known_fname), &output, false);

        assert!(result.is_err());

        // Clean up the output file
        std::fs::remove_file(Path::new(&xref_fname))?;
        std::fs::remove_file(Path::new(&known_fname))?;

        Ok(())
    }
}
