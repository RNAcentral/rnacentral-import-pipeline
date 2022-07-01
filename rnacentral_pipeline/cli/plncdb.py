# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
from furl import furl

from rnacentral_pipeline.databases.plncdb import parser
from rnacentral_pipeline.writers import entry_writer


@click.group("plncdb")
def cli():
    """
    A group of commands dealing with PLncDB data.
    """
    pass


@cli.command("parse")
@click.argument("data", type=click.Path(dir_okay=True, readable=True, file_okay=False))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def parse(data, output):
    entries = parser.parse(Path(data))
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
