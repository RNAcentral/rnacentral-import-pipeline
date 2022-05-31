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
        .with_ignore_parser_errors(true)
        .finish()
        .unwrap();

    // hstack the experiment name (derived from the filename) into the DataFrame
    let iter_exp = vec![path
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .split('-')
        .collect::<Vec<&str>>()[0..=2]
        .join("-")]
    .into_iter();
    let exp_col: Series = Series::new(
        "experiment",
        iter_exp
            .flat_map(|n| std::iter::repeat(n).take(df.height()))
            .collect::<Vec<String>>(),
    );
    // exp_col.rename("experiment");
    df.hstack_mut(&[exp_col]).unwrap();

    Ok(df)
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

        // Rename columns to remove . in the names
        let mut new_cols = Vec::new();
        for nm in input.get_column_names().iter() {
            new_cols.push(nm.replace('.', ""));
        }

        if new_cols != input.get_column_names() {
            input.set_column_names(&new_cols).unwrap_or_else(|error| {
                panic!("Failed to set column names for some reason {:?}", error);
            });
        }

        // Drop everything from the input except the genes
        let genes: DataFrame = input
            .select(select.iter())
            .unwrap_or_else(|error| match error {
                PolarsError::NotFound(_string) => {
                    panic!(
                        "{:?} was not found in the header, does the file have a header?\n{:?}",
                        select,
                        input.get_column_names()
                    );
                }
                _ => panic!("Error selecting column from input: {:?}", error),
            });

        output = match output {
            None => {
                let output_df = genes.clone();
                Some(output_df)
            }
            Some(mut output_df) => {
                // genes =/= output, so we are handling a new file
                output_df.extend(&genes).unwrap_or_else(|error| {
                    panic!("Failed to extend the output dataframe. {:?}", error)
                });
                output_df = output_df
                    .unique(None, UniqueKeepStrategy::First)
                    .unwrap_or_else(|error| {
                        panic!(
                            "Failed to parse selected column for uniqueness. {:?}",
                            error
                        )
                    });
                Some(output_df)
            }
        }
    }
    Ok(output.unwrap())
}

fn get_taxid(paths: &mut Vec<PathBuf>, lookup_table: &mut HashMap<String, String>) {
    let tax_regex = Regex::new(r".*Taxon_([0-9]{4,})").unwrap();

    while let Some(path) = paths.pop() {
        let exp_name = path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .split('-')
            .collect::<Vec<&str>>()[0..=2]
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
                        }
                        None => {}
                    };
                }
                Err(_idc) => {}
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

    let chunked_input: Vec<Vec<PathBuf>> = expr_files
        .chunks(args.chunksize)
        .map(|x| x.to_vec())
        .collect();
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
    let chunked_taxa: Vec<Vec<PathBuf>> = taxa_files
        .chunks(args.chunksize)
        .map(|x| x.to_vec())
        .collect();

    let mut lookup_table: HashMap<String, String> = HashMap::new();

    for mut tax_chunk in chunked_taxa {
        get_taxid(&mut tax_chunk, &mut lookup_table);
    }

    // Now we have the experiment - taxid lookup_table, we just need to add the column to the output df

    let expt_col = output.select(["experiment"]).unwrap();
    let mut tax_ids: Vec<String> = Vec::with_capacity(expt_col.height());
    for ex in expt_col.iter() {
        let mut warn_flag:bool = true;
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

    output
        .hstack_mut(&[Series::new("taxid", tax_ids.into_iter())])
        .unwrap();

    let out_stream: File = File::create(args.output).unwrap();
    CsvWriter::new(out_stream)
        .has_header(true)
        .finish(&mut output)
        .unwrap_or_else(|error| panic!("Something wrong writing file: {:?}", error));

    println!(
        "Processed {} in {} seconds",
        n_files,
        timer.elapsed().as_secs()
    );
}
