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

from rnacentral_pipeline.rnacentral.genes import build, write
from rnacentral_pipeline.rnacentral.genes.data import (
    Context,
    DataType,
    MemberType,
    Methods,
)
from rnacentral_pipeline.rnacentral.genes.random_forest import classify as rf_classify
from rnacentral_pipeline.rnacentral.genes.random_forest import data, preprocessing


@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


@cli.command("build")
@click.option(
    "--format",
    default="csv",
    type=click.Choice(write.Format.names(), case_sensitive=False),
)
@click.option("--method", default="rules")
@click.option("--include-genes/--no-genes", default=True)
@click.option("--include-representative/--no-representatives", default=True)
@click.option("--include-ignored/--no-ignored", default=True)
@click.option("--include-members/--no-members", default=True)
@click.option("--include-rejected/--no-rejected", default=True)
@click.option("--extended-bed/--simple-bed", default=False)
@click.argument("data_file", type=click.File("r"))
@click.argument("count_file", type=click.File("r"))
@click.argument("genes_file", type=click.Path())
@click.argument("repetitive_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def build_genes(
    data_file,
    count_file,
    genes_file,
    repetitive_file,
    output,
    format=None,
    method=None,
    include_genes=None,
    include_representative=None,
    include_members=None,
    include_ignored=None,
    include_rejected=None,
    extended_bed=None,
):
    """
    Build the genes for the given data file. The file can contain all data for a
    specific assembly.
    """

    context = Context.from_files(genes_file, repetitive_file, count_file)
    method = Methods.from_name(method)
    data = build.from_json(context, method, data_file)
    format = write.Format.from_name(format)
    allowed_members = set()
    if include_representative:
        allowed_members.add(MemberType.highlighted)
    if include_members:
        allowed_members.add(MemberType.member)

    allowed_data_types = DataType.all()
    if not include_rejected:
        allowed_data_types.remove(DataType.rejected)
    if not include_ignored:
        allowed_data_types.remove(DataType.ignored)

    output = Path(output)
    if not output.exists():
        output.mkdir()

    write.write(
        data,
        format,
        output,
        include_genes=include_genes,
        allowed_members=allowed_members,
        extended_bed=extended_bed,
        allowed_data_types=allowed_data_types,
    )


## New RF based gene models


@cli.command("fetch")
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


@cli.command("preprocess")
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
def preprocess(transcripts_file, conn_str, so_model_path, output, nearby_distance):
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
        transcripts_file, conn_str, so_model_path, nearby_distance
    )

    # Ensure output directory exists
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    features.write_parquet(output)
    click.echo(f"Generated {features.height} feature comparisons, saved to {output}")


@cli.command("classify")
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
