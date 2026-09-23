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

from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.hgnc import parser
from rnacentral_pipeline.writers import entry_writer

DATABASE = "HGNC"


@click.group("hgnc")
def cli():
    """
    Commands for HGNC data.
    """


@cli.command("map")
@click.option("--db-url", envvar="PGDATABASE")
@click.option(
    "--force-full",
    is_flag=True,
    help="Accepted for the pipeline's --force_full_import; HGNC always parses in full.",
)
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
def process_hgnc(filename, output, force_full=False, db_url=None):
    """
    Process the raw HGNC file into importable CSV files.

    Every record is mapped and written. The full signature set still goes to
    manifest.csv so a delta could start from it; see docs/incremental-parsing.md.
    """
    # HGNC is imported in full, as before delta existed: release.get_load_release_type
    # pins it to FULL, which retires every xref absent from the load, so a delta
    # parse here would retire everything unchanged. To run HGNC as a delta, lift
    # that pin and restore `manifest.load_signatures_for(db_url, DATABASE)` here.
    previous = {}
    result = parser.parse(Path(filename), db_url, previous)
    with entry_writer(Path(output)) as writer:
        writer.write(result.entries)
    manifest.write_artifacts(output, DATABASE, result.signatures, result.deletions)
