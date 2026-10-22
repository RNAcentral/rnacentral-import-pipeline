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

from rnacentral_pipeline.rnacentral.genes import convert as gff_convert
from rnacentral_pipeline.rnacentral.genes import merge as data


@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


@cli.group("utils")
def utils_cli():
    """
    Utility commands for gene processing.
    """
    pass


## Extract transcripts for the GNN gene model
@utils_cli.command("fetch")
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
    if transcripts.height > 0:
        transcripts.write_parquet(output)
    else:
        output_path.touch(exist_ok=True)
    click.echo(f"Saved {transcripts.height} transcripts to {output}")


@utils_cli.command("convert")
@click.option("--gff_file", required=True, help="Input GFF file with genes")
@click.option("--taxid", type=int, help="Taxonomy ID for the GFF file")
@click.option(
    "--regions_data",
    required=True,
    help="Path to extracted region data for this taxid, csv format",
)
def convert(gff_file, taxid, regions_data):
    """
    Convert GFF file to parquet format.

    This command reads a GFF file containing gene annotations and converts it to a
    parquet format suitable for further processing.
    """
    if not Path(gff_file).exists():
        raise click.ClickException(f"GFF file not found: {gff_file}")
    output_name = Path(gff_file).stem + "_transcripts.parquet"
    output_path = Path(gff_file).parent / output_name
    transcripts_table = gff_convert.gff_to_polars(
        Path(gff_file), taxid=taxid, regions_data=Path(regions_data)
    )
    transcripts_table.write_parquet(output_path)
    click.echo(
        f"Converted {gff_file} to {output_path} with {transcripts_table.height} transcripts."
    )


@utils_cli.command("merge")
@click.option("--previous_genes", required=True, help="Path to previous genes file")
@click.option("--next_genes", required=True, help="Path to new genes file")
@click.option("--output", required=True, help="Output file path for merged genes")
@click.option(
    "--inactive_ids",
    type=click.Path(exists=True),
    default=None,
    help="Path to inactive IDs file. Should be a CSV with one column of URS'.",
)
@click.option(
    "--prev_release_number", type=int, default=None, help="Previous release number"
)
@click.option(
    "--next_release_number", type=int, default=None, help="Next release number"
)
def merge(
    previous_genes,
    next_genes,
    output,
    inactive_ids,
    prev_release_number,
    next_release_number,
):
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
        previous_genes,
        next_genes,
        output,
        inactive_ids,
        prev_release_number,
        next_release_number,
    )

    merged_genes.write_json(output)
    click.echo(f"Merged genes saved to {output}")


@utils_cli.command("init")
@click.option("--genes", required=True, help="Path to genes file")
@click.option("--output", required=True, help="Output file path")
@click.option("--release_number", type=int, required=True, help="Release number")
def init(genes, output, release_number):
    """
    Initialize version tracking for a single release.

    Use this when only one release is available for a taxon and no merging
    is needed. Adds release tracking columns (first_release, last_release)
    and ensures consistent name formatting.
    """
    genes_path = Path(genes)
    output_path = Path(output)

    if not genes_path.exists():
        raise click.ClickException(f"Genes file not found: {genes_path}")

    initialized = data.init_genes(genes_path, release_number)
    initialized.write_json(output_path)
    click.echo(
        f"Initialized {len(initialized)} genes for release {release_number} -> {output_path}"
    )


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
@click.option(
    "--db_str",
    envvar="PGDATABASE",
)
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
@click.option(
    "--db_str",
    envvar="PGDATABASE",
)
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


@utils_cli.command("fetch-stored")
@click.option("--taxid", type=int, required=True, help="Taxonomy ID to fetch genes for")
@click.option("--output", required=True, help="Output JSON path for the stored genes")
@click.option(
    "--release_output",
    default="prev_release.txt",
    help="Path to write the previous release number (empty for a new taxon)",
)
@click.option("--db_str", envvar="PGDATABASE")
def fetch_stored(taxid, output, release_output, db_str):
    """
    Reconstruct the stored genes for a taxon from the database.

    Writes the genes as merge-format JSON (the inverse of store-genes) so they
    can be used as the previous_genes input to a forward merge. Also writes the
    previous release number to a file. An empty JSON array means no genes are
    stored yet for this taxon (a new taxon) and the caller should use init.
    """
    genes, prev_release_number = data.fetch_stored_genes(taxid, db_str)
    genes.write_json(Path(output))
    with open(release_output, "w") as handle:
        handle.write("" if prev_release_number is None else str(prev_release_number))
    click.echo(
        f"Fetched {genes.height} stored genes for taxid {taxid} "
        f"(previous release: {prev_release_number}) -> {output}"
    )


@utils_cli.command("deactivate-discarded")
@click.option("--merged", required=True, help="Path to the merged genes file")
@click.option("--taxid", type=int, required=True, help="Taxonomy ID being processed")
@click.option("--db_str", envvar="PGDATABASE")
def deactivate_discarded(merged, taxid, db_str):
    """
    Mark genes dropped by the merge as inactive in the database.

    Genes present for this taxon in the database but absent from the merged file
    have is_active set to false (their version history is preserved). Run this
    after store-genes, which re-activates the genes that are still present.
    """
    merged_path = Path(merged)
    if not merged_path.exists():
        raise click.ClickException(f"Merged genes file not found: {merged_path}")

    data.deactivate_discarded(merged_path, taxid, db_str)
    click.echo(f"Deactivated discarded genes for taxid {taxid}.")
