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

import json
from pathlib import Path

import click

from rnacentral_pipeline.databases.generic import v1
from rnacentral_pipeline.databases.mgnify.prepare import prepare_mgnify_data
from rnacentral_pipeline.writers import entry_writer


@click.group("mgnify")
def cli():
    """
    Commands for parsing MGnify data.
    """


@cli.command("parse")
@click.argument("json_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
@click.option("--db_url", envvar="PGDATABASE")
def process_json_schema(json_file, output, db_url):
    """
    Load and parse mgnify data
    """
    # MGnify json is not quite complete, so we load it up and make some modifications.
    # These include:
    # - Figuring out a taxid for each entry
    # - Adding a type based on the source model
    data = json.load(json_file)
    data = prepare_mgnify_data(data, db_url)

    entries = v1.parse(data)

    with entry_writer(Path(output)) as writer:
        writer.write(entries)
