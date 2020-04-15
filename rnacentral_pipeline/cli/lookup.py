# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.genecards_suite.core import \
    lookup as genecards

from rnacentral_pipeline.databases.intact import lookup as intact


@click.group('lookup')
def cli():
    """
    Commands for fetching data needed for genecards parsing.
    """
    pass


@cli.command('genecards')
@click.option('--db-url', envvar='PGDATABASE')
@click.argument('urs-field')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='urs-info.pickle', type=click.File('wb'))
def lookup_genecards(filename, urs_field, output, db_url):
    """
    Lookup the required information for all URS ids in the given file and write
    them to the given file.
    """
    genecards.write(filename, db_url, urs_field, output)


@cli.command('intact')
@click.option('--db-url', envvar='PGDATABASE')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='urs-info.pickle', type=click.File('wb'))
def lookup_intact(filename, output, db_url):
    intact.write(filename, db_url, output)
