use anyhow::Result;
use polars::chunked_array::builder::Utf8ChunkedBuilder;
use polars::frame::DataFrame;
use polars::prelude::IntoSeries;
use std::fs;
use std::io::Read;
use std::path::PathBuf;

use log::{info, warn};

pub fn parse_condensed_sdrf(path: &PathBuf) -> Result<DataFrame, polars::prelude::PolarsError> {
    /*
    A condensed sdrf file has 7 columns, but the last is often not delimited correctly meaning it
    is tricky to read with the polars default csv reader.

    Therefore, we will be manually parsing the file into 6 series objects (one column seems to
    always be null) and constructing a dataframe from them

    We use a chunked array builder for Utf8 strings.
    */

    info!("Loading sdrf data from {:?}", path);

    let mut file = fs::File::open(path).unwrap();
    let mut s = String::new();
    file.read_to_string(&mut s)?;

    let part_parsed: Vec<Vec<&str>> = s.lines().map(|line| line.split('\t').collect()).collect();
    let bytes_per_string: usize = 128;
    let mut exp_name = Utf8ChunkedBuilder::new("exp_name", part_parsed.len(), bytes_per_string);
    let mut assay_name = Utf8ChunkedBuilder::new("assay_name", part_parsed.len(), bytes_per_string);
    let mut feat_class = Utf8ChunkedBuilder::new("feat_class", part_parsed.len(), bytes_per_string);
    let mut feat_type = Utf8ChunkedBuilder::new("feat_type", part_parsed.len(), bytes_per_string);
    let mut feat_value = Utf8ChunkedBuilder::new("feat_value", part_parsed.len(), bytes_per_string);
    let mut ontology = Utf8ChunkedBuilder::new("ontology", part_parsed.len(), bytes_per_string);

    // There is one experiment file that does not have the empty column in line[1]
    if part_parsed.iter().map(|x| x.len()).max().unwrap() == 7 {
        for line in part_parsed.iter() {
            exp_name.append_value(line[0]);
            assay_name.append_value(line[2]); // remember line[1] will be empty
            feat_class.append_value(line[3]);
            feat_type.append_value(line[4]);
            feat_value.append_value(line[5]);
            if line.len() == 7 {
                ontology.append_value(line[6]);
            } else {
                ontology.append_null();
            }
        }
    } else {
        warn!("Unusual sdrf parsing with {} columns, not 7 for experiment {}",part_parsed.iter().map(|x| x.len()).max().unwrap(), part_parsed[0][0]);
        for line in part_parsed.iter() {
            exp_name.append_value(line[0]);
            assay_name.append_value(line[1]);
            feat_class.append_value(line[2]);
            feat_type.append_value(line[3]);
            feat_value.append_value(line[4]);
            if line.len() == 6 {
                ontology.append_value(line[5]);
            } else {
                ontology.append_null();
            }
        }
    }

    // Now have all the lines parsed with the same lengths. Try to construct a dataframe...
    DataFrame::new(vec![
        exp_name.finish().into_series(),
        assay_name.finish().into_series(),
        feat_class.finish().into_series(),
        feat_type.finish().into_series(),
        feat_value.finish().into_series(),
        ontology.finish().into_series(),
    ])
}
