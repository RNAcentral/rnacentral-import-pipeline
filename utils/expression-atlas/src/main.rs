use anyhow::Result;
use clap::Parser;
use regex::Regex;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

use polars::frame::DataFrame;
use polars::prelude::*;

use log::{info, warn};

pub mod augment;
pub mod configuration;
pub mod filtering;
pub mod sdrf;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Path where input has been copied. Must contain the config files
    #[clap(short, long, multiple_values(true))]
    input: PathBuf,

    /// An output file
    #[clap(short, long)]
    output: String,

    /// The column to select. If files have no headers, use "column_N"
    #[clap(short, long, multiple_values(true))]
    select: Vec<String>,
}

fn load_df_add_experiment(path: &PathBuf) -> Result<DataFrame, polars::prelude::PolarsError> {
    info!("Loading experiment data from {:?}", path);
    let exp_name = path
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .split('-')
        .collect::<Vec<&str>>()[0..=2]
        .join("-")
        .replace("_A", "");

    let mut exp_df: DataFrame = CsvReader::from_path(&path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .with_null_values(Some(NullValues::AllColumns("NA".to_string())))
        .infer_schema(None)
        .finish()
        .unwrap_or_else(|x| panic!("Failed on {:?} with error {:?}", path, x));

    // hstack the experiment name (derived from the filename) into the DataFrame
    let mut exp_col_arr = Utf8ChunkedBuilder::new("experiment", exp_df.height(), 128);
    for _i in 0..exp_df.height() {
        exp_col_arr.append_value(&exp_name);
    }
    // let iter_exp = std::iter::repeat([&exp_name].into_iter()).take(exp_df.height());
    let exp_col: Series = exp_col_arr.finish().into_series();
    exp_df.hstack_mut(&[exp_col]).unwrap();

    if !&exp_df.get_column_names().contains(&"GeneID") {
        // normalise column names
        info!("Standard column heading not found, normalising column names");
        info!("Column names were {:?}", &exp_df.get_column_names());
        if exp_df.get_column_names().contains(&"Gene ID") {
            exp_df.rename("Gene ID", "GeneID")?;
        } else if exp_df.get_column_names().contains(&"Gene.ID") {
            exp_df.rename("Gene.ID", "GeneID")?;
        }
        info!("Column names are now {:?}", &exp_df.get_column_names());
    }

    Ok(exp_df)
}

fn main() -> Result<()> {
    env_logger::init();
    info!("Starting Expression Atlas parser");
    let args = Args::parse();

    // Parse the config files first
    // set up the regex
    let config_re = Regex::new(r"configuration.xml").unwrap();
    let mut config_lookup: HashMap<String, configuration::Config> = HashMap::new();
    for file in fs::read_dir(&args.input)? {
        let file = file?;
        let path = file.path();
        if config_re.is_match(path.to_str().unwrap()) {
            let exp_name = path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .split('-')
                .collect::<Vec<&str>>()[0..=2]
                .join("-")
                .replace("configuration", ""); // yeah...
            let config = configuration::parse_config(&path)?;
            config_lookup.insert(exp_name, config);
        }
    }

    // Now have a hashmap with exp_name:config. We can loop over it and
    // - Check config for experiment type
    // - Construct appropriate filenames
    // - Dispatch filenames for loading, appropriate error handling if they don't exist
    // - Construct new df to merge with big one

    let differential_re = Regex::new(r".*diff.*").unwrap();

    let mut big_df = DataFrame::default();
     // String::new();
    //  data_path = PathBuf::from(&args.input);
    // let mut sdrf_path = PathBuf::from(&args.input);

    let mut gene_count: usize = 0;

    for (exp_name, config) in &config_lookup {
        let mut exp_df = DataFrame::default();
        let mut data_path = PathBuf::from(&args.input);
        let mut sdrf_path = PathBuf::from(&args.input);
        if differential_re.is_match(&config.exp_type) {
            for analysis in &config.analytics {
                let array_design = analysis.array_design.as_deref().unwrap_or("");
                if !array_design.is_empty() {
                    let data_filename: String = format!("{}_{}-analytics.tsv", exp_name, array_design);
                    data_path.push(&data_filename);

                    if !data_path.exists() {
                        warn!(
                            "File {} does not exist, skipping this experiment",
                            data_path.to_str().unwrap()
                        );
                        data_path.pop();
                        continue;
                    }

                    // Load the data
                    if exp_df.height() == 0 {
                        exp_df = load_df_add_experiment(&data_path)?;
                        data_path.pop();
                    } else {
                        exp_df = exp_df
                            .lazy()
                            .join(
                                load_df_add_experiment(&data_path)?.lazy(),
                                [col("GeneID")],
                                [col("GeneID")],
                                JoinType::Inner,
                            )
                            .select(&[col("*").exclude([
                                "Gene Name_right",
                                "experiment_right",
                                "Design Element_right",
                            ])])
                            .collect()?;
                        data_path.pop();
                    }
                } else {
                    let data_filename: String  = format!("{}-analytics.tsv", exp_name);
                    data_path.push(&data_filename);

                    if !data_path.exists() {
                        warn!(
                            "File {} does not exist, skipping this experiment",
                            data_path.to_str().unwrap()
                        );
                        data_path.pop();
                        continue;
                    }

                    exp_df = load_df_add_experiment(&data_path)?;
                    data_path.pop();
                }
            }
        } else {
            let data_filename: String = format!("{}-tpms.tsv", exp_name);
            data_path.push(&data_filename);

            if !data_path.exists() {
                warn!(
                    "File {} does not exist, skipping this experiment",
                    data_path.to_str().unwrap()
                );
                data_path.pop();
                continue;
            }

            exp_df = load_df_add_experiment(&data_path)?;
            data_path.pop();
        }

        let sdrf_filename: String = format!("{}.condensed-sdrf.tsv", exp_name);
        sdrf_path.push(&sdrf_filename);

        if !sdrf_path.exists() {
            warn!(
                "File {} does not exist, skipping this experiment",
                sdrf_path.to_str().unwrap()
            );
            sdrf_path.pop();
            data_path.pop();
            continue;
        }

        // Now load the sdrf
        let sdrf_df = sdrf::parse_condensed_sdrf(&sdrf_path)?;
        // println!("{:?}", sdrf_df);
        // filter based on differential or baseline
        if differential_re.is_match(&config.exp_type) {
            info!("Filtering experiment dataset with differential filters");
            exp_df = filtering::filter_differential(&mut exp_df)?;
            exp_df = augment::augment_differential_df(&mut exp_df, &config, &sdrf_df)?;
        } else {
            info!("Filtering with baseline filters");
            exp_df = filtering::filter_baseline(&mut exp_df);
            exp_df = augment::augment_baseline_df(&mut exp_df, &config, &sdrf_df)?;
        }

        data_path.pop();
        sdrf_path.pop();

        info!("dataframe remaining: {}", exp_df.height());
        gene_count += exp_df.height();

        // Add the newly parsed data to the big df ready for export
        if big_df.height() == 0 {
            big_df = exp_df.clone();
        } else {
            big_df.vstack_mut(&exp_df)?;
        }
    }

    info!(
        "Parsed a total of {} lines, from which {} were selected ({}%)",
        gene_count,
        big_df.height(),
        100.0 * (big_df.height() as f64) / (gene_count as f64)
    );
    println!("{:?}", big_df.height());

    info!("All files parsed, preparing to write import csvs");

    let mut output_file = fs::File::create(&args.output)?;
    JsonWriter::new(&mut output_file)
        .with_json_format(::JsonFormat::JsonLines)
        .finish(&mut big_df)?;

    Ok(())
}
