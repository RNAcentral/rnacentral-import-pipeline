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

from rnacentral_pipeline.writers import entry_writer
from rnacentral_pipeline.databases.gtrnadb import urls
from rnacentral_pipeline.databases.gtrnadb import parser


@click.group("gtrnadb")
def cli():
    """
    Commands for dealing with GtRNAdb data.
    """


@cli.command("urls-for")
@click.argument("url")
@click.argument("output", default="-", type=click.File("w"))
def urls_for(url, output):
    for url in urls.urls_for(furl(url)):
        output.write(url)
        output.write("\n")


@cli.command("parse")
@click.argument("taxonomy", type=click.Path())
@click.argument("data_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def process_gtrnadb(taxonomy, data_file, output):
    """
    Parse the GtRNAdb data with the indexed taxonomy data.
    """
    entries = parser.parse(data_file, Path(taxonomy))
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
