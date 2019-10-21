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

from rnacentral_pipeline.rnacentral import release

from rnacentral_pipeline.rnacentral import upi_ranges

from rnacentral_pipeline.databases.crs import parser as crs

from rnacentral_pipeline.rnacentral import pgloader


@click.command('run-release')
@click.option('--db_url', envvar='PGDATABASE')
def run_release(db_url=None):
    """
    A command to run the release logic in the database.
    """
    release.run(db_url)


@click.command('upi-ranges')
@click.option('--db_url', envvar='PGDATABASE')
@click.option('--table-name', default='rna')
@click.argument('chunk_size', type=int)
@click.argument('output', default='-', type=click.File('w'))
def find_upi_ranges(chunk_size, output, db_url=None, table_name=None):
    """
    This will compute the ranges to use for our each xml file in the search
    export. We want to do several chunks at once as it is faster (but not too
    man), and we want to have as large a chunk as possible. If given an a
    table_name value it will use that table, otherwise it will use the rna
    table.
    """
    upi_ranges.to_file(db_url, table_name, chunk_size, output)


@click.command('crs')
@click.argument('filename', default='-', type=click.File('r'))
@click.argument('output', default='complete_features.csv', type=click.File('w'))
def crs_data(filename, output):
    """
    This will parse the CRS file to produce a series of sequence features for
    import. The features are different from normal sequence features because
    these are 'complete', they already have a URS/taxid assigned and can just
    be inserted directly into the database.
    """
    crs.from_file(filename, output)


@click.command('validate-pgloader')
@click.argument('filename', default='-', type=click.File('r'))
def validate_pgloader(filename):
    """
    Check if pgloader ran without errors. Pgloader doesn't seem to crash when it
    should so we use this to parse the output and determine if there were any
    issues when loading. This is safer then continuing.
    """
    if not pgloader.validate(filename):
        raise click.ClickException("Pgloader produced errors")
