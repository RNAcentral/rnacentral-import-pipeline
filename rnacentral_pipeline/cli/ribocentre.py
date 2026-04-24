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

from rnacentral_pipeline.databases.generic import parser as generic
from rnacentral_pipeline.writers import entry_writer, parquet_entry_writer


@click.group("ribocentre")
def cli():
    """
    Commands for parsing ribocentre data
    """
    pass


@cli.command("parse")
@click.option("--db-url", envvar="PGDATABASE")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["csv", "parquet"]),
    default="csv",
    help="Interchange format emitted for the pgloader/DuckDB load step.",
)
@click.argument("json_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_ribocentre(json_file, output, output_format, db_url=None):
    """
    This parses the JSON that ribocentre provides us
    """
    entries = generic.parse(json_file)
    writer_factory = (
        parquet_entry_writer if output_format == "parquet" else entry_writer
    )
    with writer_factory(Path(output)) as writer:
        writer.write(entries)
