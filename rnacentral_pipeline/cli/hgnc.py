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

from rnacentral_pipeline.writers import entry_writer
from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.hgnc import parser

DATABASE = "HGNC"


@click.group("hgnc")
def cli():
    """
    Commands for HGNC data.
    """


@cli.command("map")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("filename", type=click.Path())
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_hgnc(filename, output, db_url=None):
    """
    Process the raw HGNC file into importable CSV files.

    Only records that are new or changed since the last import are mapped and
    written; dropped records go to deletions.csv and the full signature set to
    manifest.csv, both consumed by the load step. See docs/incremental-parsing.md.
    """
    previous = manifest.load_signatures_for(db_url, DATABASE)
    result = parser.parse(Path(filename), db_url, previous)
    with entry_writer(Path(output)) as writer:
        # A delta parse legitimately yields no entries when nothing changed; that is
        # the fast-path success case, not the "parser produced nothing" failure the
        # default guard catches. Safe here only because HGNC loads in DELTA mode,
        # which retires via the explicit deletions list, never by absence.
        writer.write(result.entries, allow_empty=True)
    manifest.write_artifacts(output, DATABASE, result.signatures, result.deletions)
