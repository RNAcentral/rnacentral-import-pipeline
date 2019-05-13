# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.rnacentral.text_mining import core
from rnacentral_pipeline.rnacentral.text_mining import patterns as pats


OUTPUT='mined-publications.csv'

@click.group('text-mining')
def cli():
    """
    A set of commands for extracting possible RNAcentral identifiers from text
    files. This assuems the text files are the plain text documents that
    EuropePMC provides.
    """
    pass


@cli.command('ensembl')
@click.argument('text', type=click.Path())
@click.argument('output', default=OUTPUT, type=click.File('w'))
def ensembl(text, output):
    core.write_pattern_matches(text, pats.ENSEMBL, output)


@cli.command('hgnc')
@click.argument('names', type=click.File('r'))
@click.argument('text', type=click.Path())
@click.argument('output', default=OUTPUT, type=click.File('w'))
def fixed_patterns(names, text, output):
    core.write_name_matches(text, 'hgnc', names, output)


@cli.command('mirbase')
@click.argument('text', type=click.Path())
@click.argument('output', default=OUTPUT, type=click.File('w'))
def mirbase(text, output):
    core.write_pattern_matches(text, pats.MIRBASE, output)
