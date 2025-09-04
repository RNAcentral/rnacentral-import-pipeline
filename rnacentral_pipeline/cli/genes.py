# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from pathlib import Path

import click

from rnacentral_pipeline.rnacentral.genes.random_forest import classify as rf_classify
from rnacentral_pipeline.rnacentral.genes.random_forest import convert as gff_convert
from rnacentral_pipeline.rnacentral.genes.random_forest import data, preprocessing
from rnacentral_pipeline.rnacentral.genes.random_forest import train, evaluate


@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


# Training subgroup
@cli.group("train")
def train_cli():
    """
    Commands for training gene classification models.
    """
    pass


# Inference subgroup  
@cli.group("infer")
def infer_cli():
    """
    Commands for running gene inference/classification.
    """
    pass

@cli.group("utils")
def utils_cli():
    """
    Utility commands for gene processing.
    """
    pass

## New RF based gene models

## Inference commands to fetch, preprocess and classify transcripts
@infer_cli.command("fetch")
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    required=True,
    help="Database connection string (can use PGDATABASE env var)",
)
@click.option("--taxid", type=int, required=True, help="Taxonomy ID to fetch data for")
@click.option(
    "--output",
    required=True,
    help="Output file path for transcript data (parquet format)",
)
def fetch(conn_str, taxid, output):
    """
    Fetch transcript data from the database for a given taxonomy ID.

    This command queries the database to retrieve all transcripts for the specified
    taxonomy ID and saves them to a parquet file for later processing.
    """
    click.echo(f"Fetching data for taxid {taxid}...")

    transcripts = data.fetch_data(taxid, conn_str)

    # Ensure output directory exists
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    transcripts.write_parquet(output)
    click.echo(f"Saved {transcripts.height} transcripts to {output}")


@infer_cli.command("preprocess")
@click.option(
    "--transcripts_file",
    required=True,
    help="Input parquet file containing transcript data",
)
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    help="Database connection string (required if region_ids not in transcripts)",
)
@click.option("--so_model_path", required=True, help="SO Node2Vec model path")
@click.option(
    "--output", required=True, help="Output file path for features (parquet format)"
)
@click.option(
    "--nearby_distance",
    default=1000,
    type=int,
    help="Distance threshold for identifying nearby transcripts (default: 1000)",
)
@click.option(
    "--use-parallel/--no-parallel",
    default=True,
    help="Use parallel processing across chromosome/assembly groups (default: enabled)",
)
@click.option(
    "--n-processes",
    type=int,
    help="Number of processes to use for parallel processing (default: auto-detect)",
)
def preprocess(
    transcripts_file,
    conn_str,
    so_model_path,
    output,
    nearby_distance,
    use_parallel,
    n_processes,
):
    """
    Generate features by comparing nearby transcripts.

    This command loads transcript data and generates pairwise comparison features
    for transcripts that are within the specified distance threshold.
    """
    if not Path(transcripts_file).exists():
        raise click.ClickException(f"Transcripts file not found: {transcripts_file}")

    if not Path(so_model_path).exists():
        raise click.ClickException(f"Transcripts file not found: {so_model_path}")

    click.echo(f"Loading transcripts from {transcripts_file}...")

    features = preprocessing.run_preprocessing(
        transcripts_file,
        conn_str,
        so_model_path,
        nearby_distance,
        use_parallel,
        n_processes,
    )

    # Ensure output directory exists
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    features.write_parquet(output)
    click.echo(f"Generated {features.height} feature comparisons, saved to {output}")


@infer_cli.command("classify")
@click.option(
    "--features_file", required=True, help="Input parquet file containing features"
)
@click.option(
    "--transcripts_file", required=True, help="Input parquet file containing features"
)
@click.option("--model_path", required=True, help="Path to trained model pickle file")
@click.option("--taxid", type=int, required=True, help="Taxonomy ID for gene naming")
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    required=True,
    help="Database connection string (can use PGDATABASE env var)",
)
@click.option("--output_dir", required=True, help="Output directory for results")
@click.option(
    "--seed",
    default=20240520,
    type=int,
    help="Random seed for gene naming (default: 20240520)",
)
def classify(
    features_file, transcripts_file, model_path, taxid, conn_str, output_dir, seed
):
    """
    Apply machine learning model to classify transcript pairs and generate genes.

    This command loads features, applies the trained model to classify transcript
    pairs as belonging to the same gene, then uses community detection to group
    transcripts into genes.
    """
    genes_table = rf_classify.run_final_classification(
        features_file, transcripts_file, model_path, taxid, conn_str, seed
    )
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    genes_table.write_json(output_path / f"genes_{taxid}.json")


@utils_cli.command("convert")
@click.option("--gff_file", required=True, help="Input GFF file with genes")
@click.option("--taxid", type=int, help="Taxonomy ID for the GFF file")
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    required=True,
    help="Database connection string (can use PGDATABASE env var)",
)
def convert(gff_file, taxid, conn_str):
    """
    Convert GFF file to parquet format.

    This command reads a GFF file containing gene annotations and converts it to a
    parquet format suitable for further processing.
    """
    if not Path(gff_file).exists():
        raise click.ClickException(f"GFF file not found: {gff_file}")
    output_name = Path(gff_file).stem + "_transcripts.parquet"
    output_path = Path(gff_file).parent / output_name
    transcripts_table = gff_convert.gff_to_polars(Path(gff_file), taxid=taxid, conn_str=conn_str)
    transcripts_table.write_parquet(output_path)
    click.echo(
        f"Converted {gff_file} to {output_path} with {transcripts_table.height} transcripts."
    )


@utils_cli.command("create-bedfile")
@click.option("--bed_file", required=True, help="BEDfile output path")
@click.option("--taxid", type=int, help="Taxonomy ID for the BED file")
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    required=True,
    help="Database connection string (can use PGDATABASE env var)",
)
def create_bedfile(bed_file, taxid, conn_str):
    """
    Query the genes table in the database and convert the data there into a bed file for analysis.
    """
    gff_convert.database_to_bed(
        output_path=Path(bed_file),
        taxid=taxid,
        conn_str=conn_str,
    )
    
    click.echo(
        f"Converted database records to {bed_file} with."
    )


@utils_cli.command("merge")
@click.option("--previous_genes", required=True, help="Path to previous genes file")
@click.option("--next_genes", required=True, help="Path to new genes file")
@click.option("--output", required=True, help="Output file path for merged genes")
@click.option("--inactive_ids", type=click.Path(exists=True), default=None, help="Path to inactive IDs file. Should be a CSV with one column of URS'.")
@click.option("--prev_release_number", type=int, default=None, help="Previous release number")
@click.option("--next_release_number", type=int, default=None, help="Next release number")
def merge(previous_genes, next_genes, output, inactive_ids, prev_release_number, next_release_number):
    """
    Merge two gene files, updating IDs and handling inactive genes.

    This command takes two gene files (previous and next), merges them, updates IDs,
    and handles any inactive genes based on the provided inactive IDs file.

    The merge is done based on one of two conditions:
    1. If the previous gene has a perfect overlap, the two are merged and the version is updated
    2. Within 1kb around each gene, we merge genes that have total region overlap > 0.9, and update the version

    If a gene is present in the previous release but not the next, it is discarded. Genes present in the next 
    release but not the previous are added as new genes.
    """
    previous_genes = Path(previous_genes)
    next_genes = Path(next_genes)
    output = Path(output)

    if not previous_genes.exists():
        raise click.ClickException(f"Previous genes file not found: {previous_genes}")
    if not next_genes.exists():
        raise click.ClickException(f"Next genes file not found: {next_genes}")

    merged_genes = data.merge_genes(
        previous_genes, next_genes,output, inactive_ids, prev_release_number, next_release_number
    )
    
    merged_genes.write_json(output)
    click.echo(f"Merged genes saved to {output}")


@utils_cli.command("process-metadata")
@click.argument("final_genes", type=click.Path(exists=True))
@click.argument("metadata_output", type=click.Path())
@click.option("--db_str", envvar="PGDATABASE")
def process_metadata(final_genes, metadata_output, db_str):
    metadata = data.get_metadata(final_genes, db_str)
    metadata.write_json(metadata_output)
    click.echo(f"Metadata written to {metadata_output}")




@utils_cli.command("store-genes")
@click.argument("final_genes", type=click.Path(exists=True))
@click.option("--taxid", type=int)
@click.option("--db_str", envvar="PGDATABASE",)
def store(final_genes, taxid, db_str):
    """
    Store final genes in the database.

    This command reads the final genes file and stores the gene annotations
    in the specified database.
    """
    final_genes = Path(final_genes)

    if not final_genes.exists():
        raise click.ClickException(f"Final genes file not found: {final_genes}")

    data.store_genes(final_genes, taxid, db_str)
    click.echo(f"Stored genes from {final_genes} into the database.")


@utils_cli.command("store-metadata")
@click.argument("metadata_file", type=click.Path(exists=True))
@click.option("--db_str", envvar="PGDATABASE",)
def store_metadata(metadata_file, db_str):
    """
    Store gene metadata in the database.

    This command reads the gene metadata file and stores the metadata
    in the specified database.
    """
    metadata_file = Path(metadata_file)

    if not metadata_file.exists():
        raise click.ClickException(f"Metadata file not found: {metadata_file}")

    data.store_metadata(metadata_file, db_str)
    click.echo(f"Stored metadata from {metadata_file} into the database.")


@train_cli.command("convert-csv")
@click.option("--input_csv", required=True, help="Input CSV file path")
@click.option("--output_parquet", required=True, help="Output Parquet file path")
def convert_csv(input_csv, output_parquet):
    """
    Convert the CSV dumped by psql into a parquet file for further processing
    """
    if not Path(input_csv).exists():
        raise click.ClickException(f"Input CSV file not found: {input_csv}")

    train.convert_csv_2_parquet(input_csv, output_parquet)
    click.echo(f"Converted {input_csv} to {output_parquet}")


@train_cli.command("fetch-training-data")
@click.option("--transcripts_parquet", required=True, help="The input transcripts")
@click.option("--nearby_distance", default=1000, type=int, help="Distance to search for candidate genes to compare")
@click.option("--output_parquet", required=True, help="Output Parquet file path")
def fetch_training_data(transcripts_parquet, nearby_distance, output_parquet):
    """
    Fetch training data by finding nearby transcripts and labeling them.

    This command processes the input transcripts to find nearby candidates
    within the specified distance and labels them for training.
    """
    if not Path(transcripts_parquet).exists():
        raise click.ClickException(f"Input transcripts file not found: {transcripts_parquet}")
    candidates_df = train.fetch_training_data(transcripts_parquet, nearby_distance)
    candidates_df.write_parquet(output_parquet)
    click.echo(f"Fetched training data and saved to {output_parquet}")


@train_cli.command("prepare-features")
@click.option("--candidates_parquet", required=True, help="Input Parquet file with candidates")
@click.option("--transcripts_parquet", required=True, help="Input Parquet file with transcripts")
@click.option("--so_model_path", required=True, help="SO Node2Vec model path")
@click.option("--output_parquet", required=True, help="Output Parquet file path")
def prepare_features(candidates_parquet, transcripts_parquet, so_model_path, output_parquet):
    """
    Prepare features for training the model.

    This command generates features for the candidate transcript pairs
    using the specified SO Node2Vec model.
    """
    if not Path(candidates_parquet).exists():
        raise click.ClickException(f"Input candidates file not found: {candidates_parquet}")
    if not Path(so_model_path).exists():
        raise click.ClickException(f"SO model file not found: {so_model_path}")

    features = train.build_training_features(candidates_parquet, transcripts_parquet, so_model_path)
    features.write_parquet(output_parquet)
    click.echo(f"Prepared features and saved to {output_parquet}")


@train_cli.command("split-dataset")
@click.option("--features_parquet", required=True, help="Input Parquet file with features")
@click.option("--train_path", required=True, help="Output Parquet file for training data")
@click.option("--val_path", required=True, help="Output Parquet file for validation data")
@click.option("--test_path", required=True, help="Output Parquet file for testing data")
@click.option("--test_frac", default=0.2, type=float, help="Proportion of data to use for testing (default: 0.2)")
@click.option("--val_frac", default=0.2, type=float, help="Proportion of data to use for validation (default: 0.2)")
@click.option("--seed", default=1337, type=int, help="Random seed for reproducibility (default: 1337)")
@click.option("--hub_path", default=None, help="Repo to upload dataset to on huggingface hub")
def split_dataset(features_parquet, train_path, val_path, test_path, test_frac, val_frac, seed, hub_path):
    """
    Split the dataset into training, validation, and testing sets.

    This command takes the full feature dataset and splits it into
    training, validation, and testing subsets based on the specified proportions.
    """
    if not Path(features_parquet).exists():
        raise click.ClickException(f"Input features file not found: {features_parquet}")

    train.split_datasets(features_parquet, train_path, val_path, test_path, test_frac, val_frac, seed, hub_path)
    click.echo(f"Split dataset into train ({train_path}), val ({val_path}), and test ({test_path})")


@train_cli.command("train-model")
@click.option("--train_data", required=True, help="Input Parquet file with training data")
@click.option("--model_output", required=True, help="Output folder for the trained models")
@click.option("--folds", default=5, type=int, help="Number of cross-validation folds (default: 5)")
@click.option("--exclude", multiple=True, help="Columns to exclude from features")
@click.option("--seed", default=1337, type=int, help="Random seed for reproducibility (default: 1337)")
@click.option("--basename", default="rf_model", help="Base name for output model files (default: rf_model)")
@click.option("--n_estimators", default=100, type=int, help="Number of trees in the random forest (default: 100)")
def train_model(train_data, model_output, folds, exclude, seed, basename, n_estimators):
    """
    Train a Random Forest model using the provided training data.

    This command trains a Random Forest classifier with cross-validation
    and saves the trained model to the specified output folder.
    """
    if not Path(train_data).exists():
        raise click.ClickException(f"Input training data file not found: {train_data}")
    if not Path(model_output).exists():
        Path(model_output).mkdir(parents=True, exist_ok=True)

    
    train.train(train_data, model_output, folds, exclude, seed, basename, n_estimators)
    click.echo(f"Trained models saved to {model_output}")

@train_cli.command("evaluate-model")
@click.option("--model_path", required=True, help="Path to the trained model file")
@click.option("--test_data", required=True, help="Input Parquet file with testing data")
@click.option("--metric_output", required=True, help="Output file for metrics (should be ndjson)")
@click.option("--exclude", "-e", default=None, multiple=True, help="Columns to exclude")
def evaluate_model(model_path, test_data, metric_output, exclude):
    """
    Evaluate the trained model on the test dataset.

    This command loads the trained model and evaluates its performance
    on the provided test dataset, printing out relevant metrics.
    """
    if not Path(model_path).exists():
        raise click.ClickException(f"Trained model file not found: {model_path}")
    if not Path(test_data).exists():
        raise click.ClickException(f"Input testing data file not found: {test_data}")

    excluded_columns = ["label", "comparison"]
    if exclude is not None:
        excluded_columns.extend(list(exclude))

    metrics = evaluate.evaluate_model(model_path, test_data, excluded_columns)
    click.echo("Model evaluation metrics:")
    
    click.echo(metrics)
    # Save metrics to output file
    metrics.write_ndjson(metric_output)
    
