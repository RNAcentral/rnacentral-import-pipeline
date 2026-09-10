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
import logging
import os
from pathlib import Path

import click
import psycopg2

from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ena import context, delta, parser
from rnacentral_pipeline.rnacentral.notify.slack import send_notification
from rnacentral_pipeline.writers import entry_writer

LOGGER = logging.getLogger(__name__)

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
    Emit source,accession,signature CSV for every record in an ENA .ncr chunk. Cheap:
    no ribotyper, no database. One file per chunk; the diff step unions them. The
    source label comes from the chunk's own name, which is how the diff knows which
    records a skipped source owns. See docs/incremental-parsing-ena.md.
    """
    delta.write_signatures(Path(ena_file), output)


@cli.command("source-diff")
@click.option("--db-url", envvar="PGDATABASE")
@click.option(
    "--force-full",
    is_flag=True,
    help="Ignore the stored signatures and fetch every source.",
)
@click.argument("sources_tsv", type=click.Path(exists=True))
@click.argument("to_fetch", type=click.File("w"))
@click.argument("scanned", type=click.File("w"))
@click.argument("sources_csv", type=click.File("w"))
def ena_source_diff(
    sources_tsv, to_fetch, scanned, sources_csv, force_full=False, db_url=None
):
    """
    Decide which ENA source directories are worth fetching at all.

    SOURCES_TSV is path<TAB>signature for every source in this snapshot, where the
    signature covers the names, sizes and mtimes of the source's archives. Against the
    stored signatures that gives:

      * to_fetch -- sources that are new or have changed, the only ones rsynced,
                    decompressed, split and signatured this run;
      * scanned  -- those plus sources that have gone from the snapshot: the sources
                    whose records the record-level diff may retire. Records of a
                    skipped source are not listed anywhere, and so survive untouched;
      * sources  -- database,path,signature for every current source, promoted into
                    the files table once the load has committed.

    --force-full fetches every source. It has to: a forced full run releases with
    FULL, which retires whatever is absent from the load, so skipping a source would
    retire every record in it.
    """
    current = {}
    with open(sources_tsv, "r") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            path, signature = line.split("\t")
            current[path] = signature

    stored = {} if force_full else manifest.load_file_signatures_for(db_url, DATABASE)

    changed = sorted(path for path, sig in current.items() if stored.get(path) != sig)
    vanished = sorted(path for path in stored if path not in current)

    LOGGER.info(
        "ENA sources: %d changed or new, %d gone, %d unchanged",
        len(changed),
        len(vanished),
        len(current) - len(changed),
    )

    for path in changed:
        to_fetch.write(path + "\n")
    for path in changed + vanished:
        scanned.write(path + "\n")

    sources_writer = csv.writer(sources_csv)
    for path, signature in sorted(current.items()):
        sources_writer.writerow([DATABASE, path, signature])


@cli.command("delta-diff")
@click.option("--db-url", envvar="PGDATABASE")
@click.option(
    "--force-full",
    is_flag=True,
    help="Ignore the stored manifest and parse every record.",
)
@click.argument("signatures_csv", type=click.Path(exists=True))
@click.argument("scanned_txt", type=click.Path(exists=True))
@click.argument("to_parse", type=click.File("w"))
@click.argument("deletions_csv", type=click.File("w"))
@click.argument("manifest_csv", type=click.File("w"))
def ena_delta_diff(
    signatures_csv,
    scanned_txt,
    to_parse,
    deletions_csv,
    manifest_csv,
    force_full=False,
    db_url=None,
):
    """
    Diff the collected new signatures against the stored ENA manifest and write the
    three side-channel files:

      * to_parse    -- accessions to fully parse (new or changed), or the KEEP_ALL
                       sentinel on the first delta run (no prior manifest);
      * deletions   -- database,accession rows for records that dropped out;
      * manifest    -- database,accession,signature,source for every current record.

    SIGNATURES_CSV is source,accession,signature as written per chunk, and SCANNED_TXT
    the source paths this run actually read (from `ena source-diff`). Deletion is
    restricted to those: a source skipped as unchanged contributes no signatures, so
    without the restriction every one of its records would look dropped.

    The database is only read from: the stored manifest is COPYed out and the joins
    run in polars here, so this cannot contend with anything else on the database.

    --force-full skips the diff and parses everything, for when the tracking table
    and the loaded data have drifted apart. Nothing is listed for deletion: a forced
    full run releases with FULL, which retires by absence from the load.
    """

    scanned = [line.strip() for line in open(scanned_txt) if line.strip()]
    paths = {delta.source_label(path): path for path in scanned}

    deletions_writer = csv.writer(deletions_csv)

    if force_full:
        to_parse.write(delta.KEEP_ALL + "\n")
    else:
        # Written beside the other work-directory files, never /tmp: a --contain
        # container gives /tmp a few megabytes and ENA's manifest is gigabytes.
        stored = Path("stored-manifest.csv")
        conn = psycopg2.connect(db_url)
        try:
            manifest.dump_signatures(conn, DATABASE, stored)
        finally:
            conn.close()

        result = manifest.diff_via_polars(stored, Path(signatures_csv), scanned)

        if result.is_bootstrap:
            to_parse.write(delta.KEEP_ALL + "\n")
        else:
            for accession in result.to_parse:
                to_parse.write(accession + "\n")

        for accession in result.deletions:
            deletions_writer.writerow([DATABASE, accession])

    manifest_writer = csv.writer(manifest_csv)
    with open(signatures_csv, "r", newline="") as handle:
        for label, accession, signature in csv.reader(handle):
            manifest_writer.writerow([DATABASE, accession, signature, paths[label]])


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
