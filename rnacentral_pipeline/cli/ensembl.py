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
import operator as op
from pathlib import Path

import click

from rnacentral_pipeline.databases.ensembl import parser, pseudogenes, urls
from rnacentral_pipeline.databases.ensembl.data import Division
from rnacentral_pipeline.databases.ensembl.metadata import (
    assemblies,
    compara,
    coordinate_systems,
    karyotypes,
    proteins,
)
from rnacentral_pipeline.rnacentral.notify import slack
from rnacentral_pipeline.writers import entry_writer


@click.group("ensembl")
def cli():
    """
    This is a set of commands for parsing ensembl data.
    """


@cli.command("urls-for")
@click.option("--kind", default=None)
@click.argument("division", type=click.Choice(Division.names(), case_sensitive=False))
@click.argument("ftp")
@click.argument("output", default="-", type=click.File("w"))
def vert_url(division, ftp, output, kind=None):
    """
    This is a command to generate a CSV file of urls to fetch to get Ensembl
    data from. The urls may be globs suitable for fetching with wget.
    """
    division = Division.from_name(division)
    writer = csv.writer(output, lineterminator="\n")
    rows = urls.urls_for(division, ftp)
    writer.writerows(row.writeable(kind=kind) for row in rows)


@cli.command("parse")
@click.option(
    "--family-file",
    default=None,
    type=click.Path(file_okay=True, dir_okay=False, readable=True),
)
@click.argument("division", type=click.Choice(Division.names(), case_sensitive=False))
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
def parse_data(division, embl_file, gff_file, output, family_file=None):
    """
    This will parse EMBL files from Ensembl to produce the expected CSV files.
    """
    division = Division.from_name(division)
    gff_file = Path(gff_file)
    if family_file:
        family_file = Path(family_file)
    entries = parser.parse(division, embl_file, gff_file, family_file=family_file)
    ## Send warning to slack with details about empty parse
    try:
        with entry_writer(Path(output)) as writer:
            writer.write(entries)
    except ValueError:
        print("Empty entries, implies no ncRNAs. You should check that")
        message = (f"No ncRNA entries found for {embl_file.name}, or {gff_file.name}. " 
                   + "Empty data supplied for now"
                   + ", but you should check the legitimacy of this result.\n")
        message += "For reference, the other parameters to the parser were:\n"
        message += f"division: {division}\n"
        message += f"embl_file: {embl_file.name}\n"
        message += f"gff_file: {gff_file.name}\n"
        message += f"family_file: {family_file.name}\n"

        # slack.send_notification("Ensembl parser error", message)


@cli.command("assemblies")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("connections", default="databases.json", type=click.File("r"))
@click.argument("query", default="query.sql", type=click.File("r"))
@click.argument("example_file", default="example-locations.json", type=click.File("r"))
@click.argument("known_file", default="known-assemblies.sql", type=click.File("r"))
@click.argument("output", default="assemblies.csv", type=click.File("w"))
def ensembl_write_assemblies(
    connections, query, example_file, known_file, output, db_url=None
):
    """
    This will query the ensembl databases in the connections file and write the
    output to the given file.
    """
    assemblies.write(
        connections, query, example_file, known_file, output, db_url=db_url
    )


@cli.command("coordinate-systems")
@click.argument("connections", default="databases.json", type=click.File("r"))
@click.argument("query", default="query.sql", type=click.File("r"))
@click.argument("output", default="coordinate_systems.csv", type=click.File("w"))
def ensembl_coordinates(connections, query, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    coordinate_systems.write(connections, query, output)


@cli.command("karyotypes")
@click.argument("output", default="karyotypes.csv", type=click.File("w"))
@click.argument("species", nargs=-1)
def ensembl_write_karyotypes(output, species):
    """
    Fetch all the karyotype information from all Ensembl species. This will use
    the Ensembl API to fetch the data and write to the given output file.
    """
    if not species:
        species = None
    else:
        species = set(species)
    karyotypes.write(output, species=species)


@cli.command("proteins")
@click.argument("connections", default="databases.json", type=click.File("r"))
@click.argument("query", default="query.sql", type=click.File("r"))
@click.argument("output", default="proteins.csv", type=click.File("w"))
def ensembl_proteins_cmd(connections, query, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    proteins.write(connections, query, output)


@cli.command("compara")
@click.argument("fasta", default="-", type=click.File("rb"))
@click.argument("output", default="compara.csv", type=click.File("wb"))
def ensembl_compara(fasta, output):
    """
    Parse the FASTA file of Ensembl compara data. This will produce a CSV file
    for import into the database that tracks what ensembl transcripts are
    homologous.
    """
    compara.write(fasta, output)


@cli.command("pseudogenes")
@click.argument("division", type=click.Choice(Division.names(), case_sensitive=False))
@click.argument("embl_file", type=click.File("r"))
@click.argument("output", default="ensembl-pseudogenes.csv", type=click.File("w"))
def ensembl_pseudogenes(division, embl_file, output):
    division = Division.from_name(division)
    genes = pseudogenes.parse(division, embl_file)
    genes = it.chain.from_iterable(g.writeable() for g in genes)
    writer = csv.writer(output)
    writer.writerows(genes)
