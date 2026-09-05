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

import csv
import os
from pathlib import Path

import click
import psycopg2

from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ena import context, delta, parser
from rnacentral_pipeline.rnacentral.notify.slack import send_notification
from rnacentral_pipeline.writers import entry_writer

DATABASE = "ENA"


@click.group("ena")
def cli():
    """
    Commands for parsing ENA data.
    """


@cli.command("parse")
@click.option("--counts", default="processing-results.txt")
@click.argument("ena_file", type=click.Path(file_okay=True))
@click.argument("mapping_file", type=click.Path(file_okay=True))
@click.argument("ribovore_path", type=click.Path(dir_okay=True))
@click.argument("model_lengths", type=click.Path(file_okay=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_ena(
    ena_file, mapping_file, ribovore_path, model_lengths, output, counts=None
):
    """
    Process ENA EMBL formatted files into CSV to import. The additional mapping
    file is a file containing all TPA data we are using from ENA.
    """

    ena_file = Path(ena_file)
    builder = context.ContextBuilder()
    builder.with_ribovore(Path(ribovore_path), Path(model_lengths))
    builder.with_tpa(Path(mapping_file))
    builder.with_dr(ena_file)
    with builder.context() as ctx:
        entries = parser.parse_with_context(ctx, ena_file)
        try:
            with entry_writer(Path(output)) as writer:
                writer.write(entries, allow_empty=True)
        except ValueError:
            print("No entries could be written for one of the parsed ENA files.")
            print("Sending warning to slack, but carrying on")

            # Dump this again to attach to the report
            ctx.dump_counts(Path(counts))

            message = f"No entries could be written for ENA file {ena_file}\n"
            message += "This may be correct, but you should check\n"
            message += f"Working directory: {os.getcwd()}\n"
            message += "Ribotyper log:\n"
            message += open(
                Path(ribovore_path) / "ribotyper-results.ribotyper.log", "r"
            ).read()
            message += "\n\nContext counts:\n"
            message += open(Path(counts), "r").read()

            send_notification("ENA parsing error", message)

        ctx.dump_counts(Path(counts))


@cli.command("signatures")
@click.argument("ena_file", type=click.Path(exists=True))
@click.argument("output", type=click.File("w"))
def ena_signatures(ena_file, output):
    """
    Emit accession,signature CSV for every record in an ENA .ncr chunk. Cheap: no
    ribotyper, no database. One file per chunk; the diff step unions them. See
    docs/incremental-parsing-ena.md.
    """
    delta.write_signatures(Path(ena_file), output)


@cli.command("delta-diff")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("signatures_csv", type=click.Path(exists=True))
@click.argument("to_parse", type=click.File("w"))
@click.argument("deletions_csv", type=click.File("w"))
@click.argument("manifest_csv", type=click.File("w"))
def ena_delta_diff(signatures_csv, to_parse, deletions_csv, manifest_csv, db_url=None):
    """
    Diff the collected new signatures against the stored ENA manifest, in the
    database, and write the three side-channel files:

      * to_parse    -- accessions to fully parse (new or changed), or the KEEP_ALL
                       sentinel on the first delta run (no prior manifest);
      * deletions   -- database,accession rows for records that dropped out;
      * manifest    -- database,accession,signature for every current record.
    """

    def signature_rows():
        with open(signatures_csv, "r", newline="") as handle:
            for accession, signature in csv.reader(handle):
                yield accession, signature

    conn = psycopg2.connect(db_url)
    try:
        result = manifest.diff_via_db(conn, DATABASE, signature_rows())
    finally:
        conn.close()

    if result.is_bootstrap:
        to_parse.write(delta.KEEP_ALL + "\n")
    else:
        for accession in result.to_parse:
            to_parse.write(accession + "\n")

    deletions_writer = csv.writer(deletions_csv)
    for accession in result.deletions:
        deletions_writer.writerow([DATABASE, accession])

    manifest_writer = csv.writer(manifest_csv)
    for accession, signature in signature_rows():
        manifest_writer.writerow([DATABASE, accession, signature])


@cli.command("filter")
@click.option("--only", required=True, type=click.Path(exists=True))
@click.argument("ena_file", type=click.Path(exists=True))
@click.argument("output", type=click.Path())
def ena_filter(only, ena_file, output):
    """
    Copy through only the records whose accession is listed in --only (the to-parse
    file); the KEEP_ALL sentinel copies everything. Prints the number of records
    kept so the workflow can skip ribotyper/parse for an empty result.
    """
    written = delta.filter_records(Path(ena_file), Path(only), Path(output))
    click.echo(written)
