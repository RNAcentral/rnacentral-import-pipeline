# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.genecards_suite import genecards
from rnacentral_pipeline.databases.genecards_suite import malacards
from rnacentral_pipeline.databases.genecards_suite.core import lookup

from rnacentral_pipeline.writers import entry_writer


@click.group("genecards-suite")
def cli():
    """
    A suite of commands to handle genecards
    """


@cli.command("genecards")
@click.argument("data_file", type=click.File("r"))
@click.argument("known_sequences", type=click.File("rb"))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_genecarrds(data_file, known_sequences, output):
    entries = genecards.parse(data_file, known_sequences)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)


@cli.command("malacards")
@click.argument("data_file", type=click.File("r"))
@click.argument("known_sequences", type=click.File("rb"))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_malacards(data_file, known_sequences, output):
    entries = malacards.parse(data_file, known_sequences)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)


@cli.command("lookup")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("urs-field")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="urs-info.pickle", type=click.File("wb"))
def lookup_genecards(filename, urs_field, output, db_url):
    """
    Lookup the required information for all URS ids in the given file and write
    them to the given file.
    """
    lookup.write(filename, db_url, urs_field, output)
