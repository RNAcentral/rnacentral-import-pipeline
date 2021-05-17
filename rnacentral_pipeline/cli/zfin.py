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

import json
from pathlib import Path

import click

from rnacentral_pipeline.databases import zfin
from rnacentral_pipeline.databases.generic import parser as generic
from rnacentral_pipeline.writers import entry_writer


@click.group('zfin')
def cli():
    """
    Commands for dealing with ZFIN data.
    """
    pass


@cli.command('fetch')
@click.argument('url')
@click.argument('output', default='zfin.json', type=click.File('w'))
def fetch_zfin(url, output):
    """
    Fetches ZFIN data and strips out some bad entries and fixes a few issues
    with their formatting.
    """
    json.dump(zfin.fetch(url), output)


@cli.command('parse')
@click.argument("json_file", type=click.File("r"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_json_schema(json_file, output):
    """
    This parses our JSON schema files to produce the importable CSV files.
    """
    entries = generic.parse(json_file)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
