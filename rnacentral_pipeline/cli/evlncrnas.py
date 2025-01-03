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

from rnacentral_pipeline.databases.evlncrnas.parser import parse
from rnacentral_pipeline.writers import entry_writer


@click.group("evlncrnas")
def cli():
    """
    Commands for parsing EVLncRNAs data.
    """


@cli.command("parse")
@click.argument("db_dir", type=click.Path(exists=True, dir_okay=True, readable=True))
@click.argument("db_dump", type=click.File(mode="r"), nargs=-1)
@click.option("--db-url", envvar="PGDATABASE")
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
    nargs=1,
)
def process_xlsx_files(db_dir, db_dump, output, db_url):
    """
    This parses the decompressed contents of the download into our csv files for import
    """
    entries = parse(Path(db_dir), db_dump, db_url)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
