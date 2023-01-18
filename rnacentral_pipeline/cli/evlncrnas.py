# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.evlncrnas.parser import (
    get_db_matches,
    get_ensembl_accessions,
    get_ncbi_accessions,
    parse,
    split,
)
from rnacentral_pipeline.writers import entry_writer


@click.group("evlncrnas")
def cli():
    """
    Commands for parsing EVLncRNAs data.
    """


@cli.command("download_accessions")
@click.argument(
    "accession_file", type=click.Path(exists=True, dir_okay=False, readable=True)
)
@click.argument(
    "output",
    default="ev_with_accessions.jsonl",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
)
def download_accessions(accession_file, output):
    """
    Downloads accessions from ncbi and adds them to the data from EVlncRNAs
    """
    get_accessions(Path(accession_file), Path(output))


@cli.command("download_ensembl")
@click.argument("e_file", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument(
    "output",
    default="ev_with_ensembl.jsonl",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
)
def download_accessions(e_file, output):
    """
    Downloads accessions from ncbi and adds them to the data from EVlncRNAs
    """
    get_ensembl(Path(e_file), Path(output))


@cli.command("rnc_match")
@click.argument(
    "no_accession_file", type=click.Path(exists=True, dir_okay=False, readable=True)
)
@click.argument("db_dump", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.argument(
    "output",
    default="ev_without_accessions.jsonl",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
)
def match_in_db(no_accession_file, db_dump, output):
    get_db_matches(Path(no_accession_file), Path(db_dump), Path(output))


@cli.command("split")
@click.argument("db_dir", type=click.Path(exists=True, dir_okay=True, readable=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def split_xlsx_for_accessions(db_dir, output):
    """
    Split the main excel sheet based on the presence of at least one NCBI accession for each record.

    Writes out two jsonlines files, one with and one without accessions
    """
    split(Path(db_dir), Path(output))


@cli.command("parse")
@click.argument("db_dir", type=click.Path(exists=True, dir_okay=True, readable=True))
@click.argument("db_dump", type=click.File(mode="r"))
@click.option("--db-url", envvar="PGDATABASE")
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_xlsx_files(db_dir, db_dump, output, db_url):
    """
    This parses the decompressed contents of the download into our csv files for import
    """
    entries = parse(Path(db_dir), db_dump, db_url)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
