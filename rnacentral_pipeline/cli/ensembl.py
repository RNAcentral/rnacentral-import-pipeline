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

import csv
import itertools as it
import logging
import operator as op
from pathlib import Path

import click
import psycopg2

from rnacentral_pipeline import schemas
from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ensembl import (
    delta,
    homology,
    mrna,
    parser,
    pseudogenes,
    urls,
)
from rnacentral_pipeline.databases.ensembl.metadata import (
    assemblies,
    compara,
    proteins,
)
from rnacentral_pipeline.output_format import format_option
from rnacentral_pipeline.parquet_writers import row_writer
from rnacentral_pipeline.rnacentral.notify import slack
from rnacentral_pipeline.writers import entry_writer

LOGGER = logging.getLogger(__name__)


@click.group("ensembl")
def cli():
    """
    This is a set of commands for parsing ensembl data.
    """


@cli.command("urls-for")
@click.option("--kind", default=None)
@click.option("--location", default=urls.DEFAULT_JSON_URL)
@click.argument("output", default="-", type=click.File("w"))
def vert_url(output, kind=None, location=urls.DEFAULT_JSON_URL):
    """
    Generate a CSV of (species, taxid, embl_url, gff_url) rows to fetch, one per
    species, from the unified Ensembl organisms species JSON.
    """
    writer = csv.writer(output, lineterminator="\n")
    rows = urls.urls_for(location)
    writer.writerows(row.writeable(kind=kind) for row in rows)


@cli.command("source-diff")
@click.option("--db-url", envvar="PGDATABASE")
@click.option(
    "--ftp-url", required=True, help="The FTP url the species urls start with."
)
@click.option("--nfs-path", required=True, help="Where that FTP is mirrored on /nfs.")
@click.option("--force-full", is_flag=True, help="Fetch every geneset.")
@click.option(
    "--species", multiple=True, help="Only look at these (lower-case) species."
)
@click.option("--database", default="ENSEMBL", help="Whose stored signatures to use.")
@click.argument("species_csv", type=click.Path(exists=True))
@click.argument("to_fetch", type=click.File("w"))
@click.argument("scanned", type=click.File("w"))
@click.argument("sources_csv", type=click.File("w"))
def source_diff(
    species_csv,
    to_fetch,
    scanned,
    sources_csv,
    ftp_url,
    nfs_path,
    species=(),
    force_full=False,
    db_url=None,
    database="ENSEMBL",
):
    """
    Decide which genesets are worth fetching and parsing.

    SPECIES_CSV is the `urls-for` output. A geneset whose files have the checksums it
    had at the last release is skipped, and so are its xrefs. Writes:

      * to_fetch -- the species rows that are new or changed, with their signature;
      * scanned  -- those genesets plus any that have left the list, the only ones
                    whose xrefs may be retired;
      * sources  -- database,path,signature for every current geneset, stored once
                    the release has committed.

    --force-full fetches every geneset, since a forced full release retires whatever
    the load does not contain.
    """
    with open(species_csv, newline="") as handle:
        rows = [tuple(row) for row in csv.reader(handle)]
    picked = delta.select(rows, species)
    current = delta.signatures(picked, delta.mirror(ftp_url, nfs_path))
    stored = {} if force_full else manifest.load_file_signatures_for(db_url, database)
    fetch, gone, sources = delta.plan(rows, picked, stored, current)

    LOGGER.info(
        "Ensembl: %d to fetch, %d gone, %d unchanged",
        len(fetch),
        len(gone),
        len(picked) - len(fetch),
    )
    csv.writer(to_fetch, lineterminator="\n").writerows(fetch)
    for path in [row[2] for row in fetch] + gone:
        scanned.write(path + "\n")
    csv.writer(sources_csv, lineterminator="\n").writerows(
        (database, path, sig) for path, sig in sorted(sources.items())
    )


@cli.command("delta-diff")
@click.option("--db-url", envvar="PGDATABASE")
@click.option("--database", "databases", multiple=True, default=delta.DATABASES)
@click.argument("manifest_csv", type=click.Path(exists=True))
@click.argument("scanned_txt", type=click.Path(exists=True))
@click.argument("deletions_csv", type=click.File("w"))
def delta_diff(manifest_csv, scanned_txt, deletions_csv, db_url=None, databases=()):
    """
    Write the accessions a re-parsed geneset no longer yields.

    MANIFEST_CSV is every parsed accession as database,accession,signature,source,
    and SCANNED_TXT the genesets that were read this run. Only those can lose
    xrefs; a geneset that was skipped is not in either file.
    """
    scanned = [line.strip() for line in open(scanned_txt) if line.strip()]
    conn = psycopg2.connect(db_url)
    try:
        found = delta.deletions(conn, Path(manifest_csv), scanned, Path("."), databases)
    finally:
        conn.close()
    csv.writer(deletions_csv, lineterminator="\n").writerows(found)


@cli.command("parse")
@click.option("--assembly-id", default=None, help="Use instead of the GFF's own")
@click.option("--source", default=None, help="The geneset url, to list its accessions")
@click.option("--signature", default="", help="The geneset's source signature")
@click.option(
    "--family-file",
    default=None,
    type=click.Path(file_okay=True, dir_okay=False, readable=True),
)
@click.argument("embl_file", type=click.File("r"))
@click.argument("gff_file", type=click.Path())
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
def parse_data(
    embl_file,
    gff_file,
    output,
    family_file=None,
    assembly_id=None,
    source=None,
    signature="",
):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    gff_file = Path(gff_file)
    if family_file:
        family_file = Path(family_file)
    entries = parser.parse(
        embl_file, gff_file, family_file=family_file, assembly_id=assembly_id
    )
    rows = []

    def track(entries):
        for entry in entries:
            rows.append([entry.database, entry.accession, signature, source])
            yield entry

    if source:
        entries = track(entries)
    ## Send warning to slack with details about empty parse
    try:
        with entry_writer(Path(output)) as writer:
            writer.write(entries)
    except ValueError as err:
        if str(err) not in {"Found no entries to write", "No valid entries to write"}:
            raise
        print("Empty entries, implies no ncRNAs. You should check that")
        message = (
            f"No ncRNA entries found for {embl_file.name}, or {gff_file.name}. "
            + "Empty data supplied for now"
            + ", but you should check the legitimacy of this result.\n"
        )
        message += "For reference, the other parameters to the parser were:\n"
        message += f"embl_file: {embl_file.name}\n"
        message += f"gff_file: {gff_file.name}\n"
        message += f"family_file: {family_file.name}\n"

        # slack.send_notification("Ensembl parser error", message)

    # Not *.csv: the load stream would take it for a table.
    if source:
        with (Path(output) / "accessions.manifest").open("w", newline="") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)


@cli.command("parse-mrna")
@click.option("--source", default=None, help="The geneset url, to list its accessions")
@click.option("--signature", default="", help="The geneset's source signature")
@click.argument("embl_file", type=click.File("r"))
@click.argument("gff_file", type=click.Path())
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
@format_option
def parse_mrna(embl_file, gff_file, output, source=None, signature=""):
    """
    Parse the protein coding mRNAs, and their UTRs, from an Ensembl EMBL file.
    """
    rows = []

    def track(entries):
        for entry in entries:
            rows.append((entry.database, entry.accession, signature, source))
            yield entry

    with entry_writer(Path(output)) as writer:
        writer.write(track(mrna.parse(embl_file, Path(gff_file))))
    # Not *.csv: the load stream would take it for a table.
    if source:
        with (Path(output) / "accessions.manifest").open("w", newline="") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)


@cli.command("mrna-homology")
@click.option("--homology", "homology_files", multiple=True, type=click.Path())
@click.option("--gff", "gff_files", multiple=True, type=click.Path())
@click.argument("output", default="compara.csv")
@format_option
def mrna_homology(homology_files, gff_files, output):
    """
    Link the canonical mRNAs of homologous genes, from Ensembl's homology.tsv.gz
    files, into the ensembl_compara load format.
    """
    rows = homology.rows(map(Path, homology_files), map(Path, gff_files))
    compara.write_rows(rows, output)


@cli.command("assemblies")
@click.option("--location", default=urls.DEFAULT_JSON_URL)
@click.argument("example_file", default="example-locations.json", type=click.File("r"))
@click.argument("ucsc_file", default="ucsc.json", type=click.File("r"))
@click.argument("output", default="assemblies.csv", type=click.Path())
@format_option
def ensembl_write_assemblies(example_file, ucsc_file, output, location):
    """
    Write the ensembl_assembly rows for the genesets selected from species.json.
    """
    assemblies.write(location, example_file, ucsc_file, Path(output))


@cli.command("proteins")
@click.argument("connections", default="databases.json", type=click.File("r"))
@click.argument("query", default="query.sql", type=click.File("r"))
@click.argument("output", default="proteins.csv", type=click.Path())
@format_option
def ensembl_proteins_cmd(connections, query, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    proteins.write(connections, query, Path(output))


@cli.command("compara")
@click.argument("fasta", default="-", type=click.File("rb"))
@click.argument("output", default="compara.csv")
@format_option
def ensembl_compara(fasta, output):
    """
    Parse the FASTA file of Ensembl compara data. This will produce a CSV file
    for import into the database that tracks what ensembl transcripts are
    homologous.
    """
    compara.write(fasta, output)


@cli.command("pseudogenes")
@click.argument("embl_file", type=click.File("r"))
@click.argument("output", default="ensembl-pseudogenes.csv", type=click.Path())
@format_option
def ensembl_pseudogenes(embl_file, output):
    genes = pseudogenes.parse(embl_file)
    genes = it.chain.from_iterable(g.writeable() for g in genes)
    with row_writer(Path(output), schemas.ENSEMBL_PSEUDOGENES) as writer:
        writer.writerows(genes)
