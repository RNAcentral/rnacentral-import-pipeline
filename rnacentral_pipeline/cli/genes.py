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

from rnacentral_pipeline.rnacentral.genes import split
from rnacentral_pipeline.rnacentral.genes import build

@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


@cli.command('split-coordinates')
@click.argument('data_file', type=click.File('r'))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def split_coordinates(data_file, output):
    """
    Take a json file that represents all sequences to possibly into genes in a
    species and split the data by RNA type and chromosomes. Each individual file
    can then be processed separately.
    """
    split.write_split(data_file, output)


@cli.command('build')
@click.argument('data_file', type=click.File('r'))
@click.argument('output', type=click.File('2'))
def build(data_file, output):
    """
    Build the genes for the given data file. This assumes that the file is
    already split into reasonable chunks.
    """
    build.write_genes(data, output)
