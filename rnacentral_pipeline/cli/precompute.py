# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline import writers
from rnacentral_pipeline.output_format import format_option, is_parquet
from rnacentral_pipeline.rnacentral.precompute import extract as pre_extract
from rnacentral_pipeline.rnacentral.precompute import metadata as pre_metadata
from rnacentral_pipeline.rnacentral.precompute import normalize as pre_normalize
from rnacentral_pipeline.rnacentral.precompute import process as pre
from rnacentral_pipeline.rnacentral.precompute import ranges as pre_ranges
from rnacentral_pipeline.rnacentral.precompute import select_releases as pre_select


@click.group("precompute")
def cli():
    """
    This is a group of commands for dealing with our precompute steps.
    """
    pass


@cli.command("from-file")
@click.argument("context", type=click.Path(dir_okay=True, file_okay=True))
@click.argument("json_file", type=click.Path(dir_okay=False, file_okay=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@format_option
def precompute_from_file(context, json_file, output):
    """
    This command will take the output produced by the precompute query and
    process the results into a CSV that can be loaded into the database.
    """
    updates = pre.parse(Path(context), Path(json_file))
    out_path = Path(output)
    if is_parquet():
        opener = pre.parquet_writer(out_path)
    else:
        opener = writers.build(pre.Writer, out_path)
    with opener as writer:
        writer.write(updates)


@cli.command("normalize")
@click.argument("accessions", type=click.Path(dir_okay=False, file_okay=True))
@click.argument("metadata", type=click.Path(dir_okay=False, file_okay=True))
@click.argument(
    "output", type=click.Path(writable=True, dir_okay=False, file_okay=True)
)
def precompute_normalize(accessions, metadata, output):
    """
    Join a raw, ungrouped accessions file (typed Parquet produced by
    `rnac precompute extract-query` from get-accessions/query.sql) with a
    merged metadata.json (now produced by `rnac precompute metadata-build`),
    producing the same merged/normalized output the Rust `precompute
    normalize` binary does. Replaces that binary and the separate
    `precompute group-accessions` step (grouping now happens here).
    """
    pre_normalize.write(Path(accessions), Path(metadata), Path(output))


@cli.command("metadata-build")
@click.argument("ranges", type=click.Path(dir_okay=False, file_okay=True))
@click.argument("raw_dir", type=click.Path(dir_okay=True, file_okay=False))
@click.argument(
    "output_dir", type=click.Path(writable=True, dir_okay=True, file_okay=False)
)
def precompute_metadata_build(ranges, raw_dir, output_dir):
    """
    Replace the Rust `precompute metadata group`+`merge` binaries. Reads
    the raw per-type metadata files in RAW_DIR (basic.parquet,
    coordinates.parquet, rfam-hits.parquet, r2dt-hits.parquet,
    previous.parquet, orfs.parquet, stopfree.parquet, tcode.parquet -
    produced by `rnac precompute extract-query`) and RANGES (urs_taxid.csv,
    upi_min/ut_min/ut_max rows), and writes one range-scoped metadata file
    per row into OUTPUT_DIR, plus a manifest.csv mapping upi_min to each
    file's path.
    """
    pre_metadata.build(Path(ranges), Path(raw_dir), Path(output_dir))


@cli.command("extract-query")
@click.argument("sql_file", type=click.Path(dir_okay=False, file_okay=True))
@click.argument(
    "output", type=click.Path(writable=True, dir_okay=False, file_okay=True)
)
@click.option("--range", "range_", type=(int, int), default=None)
@click.option("--pg-url", envvar="PGDATABASE")
def precompute_extract_query(sql_file, output, range_, pg_url):
    """
    Run a query from files/precompute/queries/*.sql (or
    get-accessions/query.sql) against Postgres via connectorx and write the
    typed result to Parquet. Replaces piping psql's json_build_object(...)
    output to a file. --range substitutes real bounds for the query's own
    `:min`/`:max` placeholders; omitting it substitutes NULL for both,
    matching the query's own "unranged" case.
    """
    pre_extract.extract_query(Path(sql_file), Path(output), pg_url, range=range_)


@cli.command("select-outdated")
@click.argument("xref", type=click.Path(dir_okay=False, file_okay=True))
@click.argument("known", type=click.Path(dir_okay=False, file_okay=True))
@click.argument(
    "output", type=click.Path(writable=True, dir_okay=False, file_okay=True)
)
def precompute_select_outdated(xref, known, output):
    """
    Select which urs_taxid pairs need (re)computing, by comparing each
    pair's xref release against its last-known-precomputed release (typed
    Parquet, both produced by `rnac precompute extract-query`). Replaces
    the Rust `precompute select` binary and the external `sort` steps its
    sorted-merge-join needed.
    """
    pre_select.write(Path(xref), Path(known), Path(output))


@cli.command("upi-taxid-ranges")
@click.option("--db-url", envvar="PGDATABASE")
@click.option("--tablename", default="precompute_urs_taxid")
@click.argument("ranges", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def precompute_find_upi_ranges(ranges, output, db_url=None, tablename=None):
    pre_ranges.write(ranges, output, db_url=db_url, tablename=tablename)
