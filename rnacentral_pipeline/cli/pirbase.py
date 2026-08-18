# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.pirbase import data as pirbase_data
from rnacentral_pipeline.databases.pirbase import parser
from rnacentral_pipeline.writers import entry_writer


@click.group("pirbase")
def cli():
    """
    A group of commands dealing with piRBase data.
    """
    pass


@cli.command("urls")
@click.argument("base")
@click.argument("output", default="-", type=click.File("w"))
def urls(base, output):
    """
    Write every file to fetch as code,kind,url. Driven by the species table so
    the download list and the taxids cannot drift apart.
    """

    for code, kind, url in pirbase_data.urls(base.rstrip("/")):
        output.write(f"{code},{kind},{url}\n")


@cli.command("build-known")
@click.argument("md5s", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(writable=True, dir_okay=False))
def build_known(md5s, output):
    """
    Turn the md5 dump into the sqlite index the full sets are filtered against.
    Built once per run; every species then opens it read only.
    """

    parser.build_known(Path(md5s), Path(output))


@cli.command("parse")
@click.argument("code")
@click.argument("fasta", type=click.Path(exists=True, dir_okay=False))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
@click.option(
    "--known",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="md5s already in RNAcentral; omit to keep every sequence, as the gold sets need.",
)
def parse(code, fasta, output, known):
    species = pirbase_data.species(code)
    with open(fasta, "r") as handle:
        entries = parser.parse(
            handle, species, known_path=Path(known) if known else None
        )
        with entry_writer(Path(output)) as writer:
            writer.write(entries)
