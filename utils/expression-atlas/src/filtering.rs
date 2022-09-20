/// This module implements the filtering that EA uses in their webapp
/// NB: This is all reverse engineered
use anyhow::Result;
use polars::frame::DataFrame;
use polars::prelude::*;
use polars::series::Series;
use regex::Regex;

fn lowercase_fn(ent: &Option<&str>) -> Option<String> {
    ent.as_ref().map(|ent| ent.to_lowercase())
    // match ent {
    //     None => None,
    //     Some(ent) => Some(ent.to_lowercase()),
    // }
}

fn fix_bad_infinities(str_val: &Series) -> Series {
    let lowercased = str_val
        .utf8()
        .unwrap()
        .into_iter()
        .map(|x| lowercase_fn(&x))
        .collect::<Vec<Option<String>>>()
        .into_iter()
        .map(|x| x.map(|x| x.parse::<f64>().unwrap()))
        .collect::<Vec<Option<f64>>>();

    Series::from_iter(lowercased)
}

fn baseline_get_median_gt_zero(str_val: &Series) -> Series {
    let lists = str_val.utf8().unwrap().into_iter().map(|x| {
        x.unwrap().split(',').into_iter().map(|y| y.parse::<f64>().unwrap()).collect::<Vec<f64>>()
    });

    let medians: Vec<bool> =
        lists.into_iter().map(|x| Series::from_iter(x).median().unwrap() > 0.0).collect();
    Series::from_iter(medians)
}

/// This function will filter the differential results based on:
/// - non-null p value
/// - absolute log2 fold change greater than 1
// find the p value and log fold columns
pub fn filter_differential(input: &DataFrame) -> Result<DataFrame> {
    let pv_regex = Regex::new(r".*p-value.*").unwrap();
    let log_fold_regex = Regex::new(r".*log2.*").unwrap();

    let mut inter = input.clone();

    for column in input.get_column_names_owned() {
        // Check for badly parsed infinities
        if pv_regex.is_match(&column) {
            if !inter.column(&column)?.dtype().is_numeric() {
                inter.apply(&column, fix_bad_infinities)?;
            }
            inter = inter
                .lazy()
                .with_column(
                    when(col(&column).lt(lit(0.05f64)))
                        .then(lit(true).alias(&column))
                        .otherwise(lit(false).alias(&column)),
                )
                .collect()?;
        } else if log_fold_regex.is_match(&column) {
            if !inter.column(&column)?.dtype().is_numeric() {
                inter.apply(&column, fix_bad_infinities)?;
            }
            inter = inter
                .lazy()
                .with_column(
                    when(col(&column).abs().gt_eq(lit(1.0f64)))
                        .then(lit(true).alias(&column))
                        .otherwise(lit(false).alias(&column)),
                )
                .collect()?;
        }
    }
    Ok(inter)
}

pub fn filter_baseline(input: &mut DataFrame) -> DataFrame {
    // This is a baseline experiment.
    // Find columns starting with lower case g, then apply the function to convert to
    // medain and select greater than zero
    let mut meas = Vec::<Expr>::new();

    for column in input.get_column_names_owned() {
        if column.starts_with('g') {
            input.apply(&column, baseline_get_median_gt_zero).unwrap();
            meas.push(col(&column));
        }
    }
    // Selection should now have all the gN column names in it
    input
        .clone()
        .lazy()
        .filter(
            any_exprs(&meas[0..meas.len() / 2]).or(any_exprs(&meas[meas.len() / 2..meas.len()])),
        ) // If we try to do the whole thing at once, we get a stack overflow
        .collect()
        .unwrap()
}
