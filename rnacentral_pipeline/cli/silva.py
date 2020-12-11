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

import click

from rnacentral_pipeline.databases.ncbi import taxonomy
from rnacentral_pipeline.databases.silva import parser as silva


@click.group("silva")
def cli():
    """
    Commands for dealing with SILVA data.
    """


@cli.command("parse")
@click.argument("silva-file", type=click.File("r"))
@click.argument("taxonomy", type=click.Path())
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False,),
)
def process_silva(silva_file, taxonomy output):
    write_entries(silva.parse, output, silva_file, taxonomy)



@cli.command('index-taxonomy')
@click.argument('ncbi', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
@click.argument('output', type=click.Path())
def index_taxonomy(ncbi, output):
    taxonomy.index(ncbi, output)
