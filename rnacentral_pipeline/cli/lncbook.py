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

from rnacentral_pipeline.databases.lncbook import parser
from rnacentral_pipeline.writers import entry_writer


@click.group('lncbook')
def cli():
    """
    Commands for fetching data needed for genecards parsing.
    """
    pass


@cli.command("parse")
@click.argument("json_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_lncbook(json_file, output):
    """
    This parses our JSON files from LncBook to produce the importable CSV files.
    LncBook has some special logic around sequences (only include hg38), that
    other JSON databases do not have.
    """
    entries = parser.parse(json_file)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
