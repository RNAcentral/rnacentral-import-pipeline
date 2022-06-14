use anyhow::Result;
use clap::Parser;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::PathBuf;
use std::time::Instant;

use polars::frame::DataFrame;
use polars::prelude::*;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// List of input files. Ideally with headers
    #[clap(short, long, multiple_values(true))]
    input: Vec<PathBuf>,

    /// An output file
    #[clap(short, long)]
    output: String,

    /// The column to select. If files have no headers, use "column_N"
    #[clap(short, long, multiple_values(true))]
    select: Vec<String>,

    /// The size of chunks to process at once
    #[clap(short, long, default_value_t = 100)]
    chunksize: usize,
}

fn load_data(path: &PathBuf, delim: u8) -> Result<DataFrame> {
    let mut df: DataFrame = CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(delim)
        .with_null_values(Some(NullValues::AllColumns("NA".to_string())))
        .infer_schema(None)
        .finish()
        .unwrap();

    // hstack the experiment name (derived from the filename) into the DataFrame
    let iter_exp =
        vec![path.file_name().unwrap().to_str().unwrap().split('-').collect::<Vec<&str>>()[0..=2]
            .join("-")]
        .into_iter();
    let exp_col: Series = Series::new(
        "experiment",
        iter_exp.flat_map(|n| std::iter::repeat(n).take(df.height())).collect::<Vec<String>>(),
    );
    // exp_col.rename("experiment");
    df.hstack_mut(&[exp_col]).unwrap();

    Ok(df)
}

fn baseline_get_median_gt_zero(str_val: &Series) -> Series {
    let lists = str_val.utf8().unwrap().into_iter().map(|x| {
        x.unwrap().split(',').into_iter().map(|y| y.parse::<f64>().unwrap()).collect::<Vec<f64>>()
    });

    let medians: Vec<bool> =
        lists.into_iter().map(|x| Series::from_iter(x).median().unwrap() > 0.0).collect();
    Series::from_iter(medians)
}

fn filter_input(input: &mut DataFrame, baseline: bool) -> DataFrame {
    if baseline {
        filter_baseline(input)
    } else {
        filter_differential(input)
    }
}

fn filter_baseline(input: &mut DataFrame) -> DataFrame {
    // This is a baseline experiment.
    // Find columns starting with lower case g, then apply the function to convert to
    // medain and select greater than zero
    let mut meas = Vec::<String>::new();

    for col in input.get_column_names_owned() {
        if col.starts_with('g') {
            input.apply(&col, baseline_get_median_gt_zero).unwrap();
            meas.push(col);
        }
    }
    println!("{:?}", meas);
    // Selection should now have all the gN column names in it
    input
        .clone()
        .lazy()
        .filter(any_exprs(meas.into_iter().map(|x| col(&x)).collect::<Vec<Expr>>()))
        .collect()
        .unwrap()
}

fn filter_differential(input: &DataFrame) -> DataFrame {
    // find the p value and log fold columns
    let pv_regex = Regex::new(r".*p-value.*").unwrap();
    let log_fold_regex = Regex::new(r".*log2.*").unwrap();
    let mut exprs: Vec<Expr> = Vec::new();
    let mut pv: Expr;
    let mut lf: Expr;
    for column in input.get_column_names_owned() {
        if pv_regex.is_match(&column) {
            pv = col(&column).is_not_null();
            exprs.push(pv);
        } else if log_fold_regex.is_match(&column) {
            lf = col(&column).neq(lit(0));
            exprs.push(lf);
        }
    }

    input.clone().lazy().filter(all_exprs(exprs)).collect().unwrap()
}

fn load_chunk(
    paths: &mut Vec<PathBuf>,
    select: &mut Vec<String>,
) -> Result<DataFrame, anyhow::Error> {
    /*
    Work with the chunks of input to reduce them into a single dataframe of the genes we want from
    that particular set of experiments. Should return a single dataframe
    */
    let mut output: Option<DataFrame> = None;
    while let Some(infile) = paths.pop() {
        // Load the input Csv
        let mut input: DataFrame = if infile.extension().unwrap() == "tsv" {
            load_data(&infile, b'\t')
        } else {
            load_data(&infile, b',')
        }
        .unwrap_or_else(|error| match error.downcast_ref::<PolarsError>() {
            Some(PolarsError::Io(_string)) => panic!("Input file does not exist! {:?}", infile),
            _ => panic!("An error occurred! {:?}", error),
        });

        // Rename columns to remove . in the names. Now also remove spaces
        let mut new_cols = Vec::new();
        let mut c: u8 = 0;
        for nm in input.get_column_names().iter() {
            let mut new_name = nm.replace('.', "").replace(' ', "");
            if new_cols.contains(&new_name) {
                c += 1;
                new_name += &c.to_string(); // Hopefully avoid duplicate names
            }
            new_cols.push(new_name);
        }

        if new_cols != input.get_column_names() {
            input.set_column_names(&new_cols).unwrap_or_else(|error| {
                panic!("Failed to set column names for some reason {:?}", error);
            });
        }

        // Need to detect which type of experiment it is - differential or baseline
        let type_re = Regex::new(r"tpms").unwrap();
        let filtered = filter_input(&mut input, type_re.is_match(infile.to_str().unwrap()));
        // Drop everything from the input except the genes
        let genes: DataFrame = filtered.select(select.iter()).unwrap_or_else(|error| match error {
            PolarsError::NotFound(_string) => {
                panic!(
                    "{:?} was not found in the header, does the file have a header?\n{:?}",
                    select,
                    input.get_column_names()
                );
            },
            _ => panic!("Error selecting column from input: {:?}", error),
        });

        output = match output {
            None => {
                let output_df = genes.clone();
                Some(output_df)
            },
            Some(mut output_df) => {
                // genes =/= output, so we are handling a new file
                output_df.extend(&genes).unwrap_or_else(|error| {
                    panic!("Failed to extend the output dataframe. {:?}", error)
                });
                output_df =
                    output_df.unique(None, UniqueKeepStrategy::First).unwrap_or_else(|error| {
                        panic!("Failed to parse selected column for uniqueness. {:?}", error)
                    });
                Some(output_df)
            },
        }
    }
    Ok(output.unwrap())
}

fn get_taxid(paths: &mut Vec<PathBuf>, lookup_table: &mut HashMap<String, String>) {
    let tax_regex = Regex::new(r".*Taxon_([0-9]{4,})").unwrap();

    while let Some(path) = paths.pop() {
        let exp_name =
            path.file_name().unwrap().to_str().unwrap().split('-').collect::<Vec<&str>>()[0..=2]
                .join("-")
                .replace(".condensed", "");

        let file = File::open(path).unwrap();
        let lines = io::BufReader::new(file).lines();
        for line in lines {
            match line {
                Ok(line) => {
                    let captures = tax_regex.captures(line.as_str());
                    match captures {
                        Some(caps) => {
                            let taxid = caps.get(1).map_or("", |m| m.as_str());
                            lookup_table.insert(exp_name, taxid.to_string());
                            break;
                        },
                        None => {},
                    };
                },
                Err(_idc) => {},
            }
        }
    }
}

fn main() {
    let timer = Instant::now();

    let mut args = Args::parse();

    // Always select the experiment
    args.select.push("experiment".to_string());

    // separate expression data from taxa data
    let mut taxa_files: Vec<PathBuf> = Vec::new();
    let mut expr_files: Vec<PathBuf> = Vec::new();

    let re = Regex::new(r"sdrf.tsv").unwrap();

    for infile in args.input {
        if re.is_match(infile.to_str().unwrap()) {
            // File is experiment setup data, parse for taxa
            taxa_files.push(infile);
        } else {
            // Should only be expression data, parse for genes
            expr_files.push(infile);
        }
    }

    println!("{} taxa files to parse", taxa_files.len());
    println!("{} expression files to parse", expr_files.len());

    let n_files = expr_files.len();

    let mut dataframes: Vec<DataFrame> = Vec::new();

    let chunked_input: Vec<Vec<PathBuf>> =
        expr_files.chunks(args.chunksize).map(|x| x.to_vec()).collect();
    // Read everything into a big vector
    let mut n_chunks = 0;
    for mut files_chunk in chunked_input {
        dataframes.push(
            load_chunk(&mut files_chunk, &mut args.select)
                .unwrap_or_else(|_error| panic!("Something wrong in one of the reads, aborting")),
        );
        n_chunks += 1;
        println!("Done {} files", n_chunks * args.chunksize);
    }

    let mut output: DataFrame = dataframes[0].clone();

    for item in dataframes.iter() {
        output.extend(item).unwrap_or_else(|error| {
            panic!("Unable to extend the DataFrame! Out of memory? {:?}", error)
        });
        output.unique(None, UniqueKeepStrategy::First).unwrap();
    }

    // Now parse the experiment files to get taxa
    let chunked_taxa: Vec<Vec<PathBuf>> =
        taxa_files.chunks(args.chunksize).map(|x| x.to_vec()).collect();

    let mut lookup_table: HashMap<String, String> = HashMap::new();

    for mut tax_chunk in chunked_taxa {
        get_taxid(&mut tax_chunk, &mut lookup_table);
    }

    // Now we have the experiment - taxid lookup_table, we just need to add the column to the output df

    let expt_col = output.select(["experiment"]).unwrap();
    let mut tax_ids: Vec<String> = Vec::with_capacity(expt_col.height());
    for ex in expt_col.iter() {
        let mut warn_flag: bool = true;
        let uhy = ex.utf8().unwrap();
        for (idx, x) in uhy.into_iter().enumerate() {
            if lookup_table.contains_key(x.unwrap()) {
                tax_ids.insert(idx, lookup_table.get(x.unwrap()).unwrap().to_string());
            } else {
                tax_ids.insert(idx, String::new());
                if warn_flag {
                    println!("Experiment {} does not name a taxon", x.unwrap());
                    warn_flag = false;
                }
            }
        }
    }

    output.hstack_mut(&[Series::new("taxid", tax_ids.into_iter())]).unwrap();

    let out_stream: File = File::create(args.output).unwrap();
    CsvWriter::new(out_stream)
        .has_header(true)
        .finish(&mut output)
        .unwrap_or_else(|error| panic!("Something wrong writing file: {:?}", error));

    println!("Processed {} in {} seconds", n_files, timer.elapsed().as_secs());
}
