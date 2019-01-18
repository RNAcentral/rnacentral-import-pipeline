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

from rnacentral_pipeline.databases.ensembl.metadata import proteins
from rnacentral_pipeline.databases.ensembl.metadata import assemblies
from rnacentral_pipeline.databases.ensembl.metadata import karyotypes
from rnacentral_pipeline.databases.ensembl.metadata import coordinate_systems


@click.group('ensembl')
def cli():
    """
    This is a set of commands for dealing with processing protein information.
    We don't have much in the way of protein summary but sometimes we do need a
    little for display.
    """
    pass


@cli.command('assemblies')
@click.option('--db_url', envvar='PGDATABASE')
@click.argument('connections', default='databases.json', type=click.File('rb'))
@click.argument('query', default='query.sql', type=click.File('rb'))
@click.argument('example_file', default='example-locations.json', type=click.File('rb'))
@click.argument('example_query', default='find-examples.sql', type=click.File('r'))
@click.argument('output', default='assemblies.csv', type=click.File('wb'))
def ensembl_write_assemblies(connections, query, example_file, example_query, output,
                             db_url=None):
    """
    This will query the ensembl databases in the connections file and write the
    output to the given file.
    """
    assemblies.write(connections, query, example_file, example_query, db_url, output)


@cli.command('coordinate-systems')
@click.argument('connections', default='databases.json', type=click.File('rb'))
@click.argument('query', default='query.sql', type=click.File('rb'))
@click.argument('output', default='coordinate_systems.csv', type=click.File('wb'))
def ensembl_coordinates(connections, query, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    coordinate_systems.write(connections, query, output)


@cli.command('karyotypes')
@click.argument('output', default='karyotypes.csv', type=click.File('wb'))
def ensembl_write_karyotypes(output):
    """
    Fetch all the karyotype information from all Ensembl species. This will use
    the Ensembl API to fetch the data and write to the given output file.
    """
    karyotypes.write(output)


@cli.command('proteins')
@click.argument('connections', default='databases.json', type=click.File('rb'))
@click.argument('query', default='query.sql-', type=click.File('rb'))
@click.argument('output', default='proteins.csv', type=click.File('wb'))
def ensembl_proteins_cmd(connections, query, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    proteins.write(connections, query, output)
