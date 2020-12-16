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
from pathlib import Path

from rnacentral_pipeline.databases.ncbi import taxonomy


@click.group("context")
def cli():
    """
    Commands for building the context needed for parsing data.
    """


@cli.command('build')
@click.argument('ncbi', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
@click.argument('output', type=click.Path())
def index_taxonomy(ncbi, output):
    taxonomy.index(Path(ncbi), output)
