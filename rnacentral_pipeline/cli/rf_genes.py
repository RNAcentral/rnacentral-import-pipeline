# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.random_forest_genes import (
    convert,
    extract,
    preprocess,
    split,
    train,
)


@click.group("rf-genes")
def cli():
    """
    Commands to train and apply random forest transcript->gene models
    """
    pass


@cli.command("convert")
@click.argument("csv_file", type=click.File("r"))
@click.argument("parquet_file", type=click.File("w"))
def convert_csv_2_parquet(csv_file, parquet_file):
    convert.csv_2_parquet(csv_file.name, parquet_file.name)


@cli.command("extract")
@click.argument("gene_dataframe", type=click.File("r"))
@click.argument("taxid", type=int)
@click.option(
    "--output_transcripts", help="Where to write the resulting dataframe", default=None
)
@click.option(
    "--output_genes", help="Where to write the resulting dataframe", default=None
)
@click.option("--ensembl_only", is_flag=True, default=False)
@click.option(
    "--pre_gene_set", default=None, help="A predownloade ensembl gene + coordinate set"
)
@click.option("--nearby_distance", default=1e3)
def extract_candidates_transcripts(
    gene_dataframe,
    taxid,
    output_transcripts,
    output_genes,
    ensembl_only,
    pre_gene_set,
    nearby_distance,
):
    extract.genes_transcripts(
        gene_dataframe,
        taxid,
        output_transcripts,
        output_genes,
        ensembl_only,
        pre_gene_set,
        nearby_distance,
    )


@cli.command("preprocess")
@click.argument("candidates", type=click.File("r"))
@click.argument("transcript_data", type=click.File("r"))
@click.argument("feature_output", type=click.File("w"))
def preprocess_transcripts(candidates, transcript_data, feature_output):
    preprocess.transcripts(candidates.name, transcript_data.name, feature_output.name)


@cli.command("split")
@click.argument("input_data")
@click.argument("train_path")
@click.argument("val_path")
@click.argument("test_path")
@click.option("--test_frac", type=float, default=0.2)
@click.option("--val_frac", type=float, default=0.2)
@click.option("--seed", type=int, default=1337)
@click.option("--hub_repo", type=str, default=None)
def split_dataset(
    input_data, train_path, test_path, val_path, test_frac, val_frac, seed, hub_repo
):
    split.train_test_val(
        input_data, train_path, test_path, val_path, test_frac, val_frac, seed, hub_repo
    )


@cli.command("train_nokf")
@click.argument("training_data")
@click.argument("model_file")
@click.option("--n_estimators", type=int, default=100)
@click.option("--exclude", "-e", type=str, multiple=True)
@click.option("--seed", "-s", type=int, default=42)
@click.option("--hub_path", type=str, default=None)
def train_no_kfold(training_data, model_file, n_estimators, exclude, seed, hub_path):
    train.no_kfold(training_data, model_file, n_estimators, exclude, seed, hub_path)


@cli.command("train_kf")
@click.argument("training_data")
@click.argument("output_folder")
@click.option("--folds", type=int, default=5, help="How many folds to run")
@click.option("--exclude", "-e", multiple=True, help="Columns to exclude")
@click.option("--seed", type=int, default=1337, help="Random number seed")
@click.option("--basename", default="model_fold_", help="Basename for the model file")
@click.option("--n_estimators", type=int, default=100)
def train_with_kfold(
    training_data, output_folder, folds, exclude, seed, basename, n_estimators
):
    train.with_kfold(
        training_data, output_folder, folds, exclude, seed, basename, n_estimators
    )
