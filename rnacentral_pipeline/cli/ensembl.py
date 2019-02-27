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

from rnacentral_pipeline.databases.ensembl.metadata import assemblies
from rnacentral_pipeline.databases.ensembl.metadata import compara
from rnacentral_pipeline.databases.ensembl.metadata import coordinate_systems
from rnacentral_pipeline.databases.ensembl.metadata import karyotypes
from rnacentral_pipeline.databases.ensembl.metadata import proteins


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
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql', type=click.File('r'))
@click.argument('example_file', default='example-locations.json', type=click.File('r'))
@click.argument('example_query', default='find-examples.sql', type=click.File('r'))
@click.argument('output', default='assemblies.csv', type=click.File('w'))
def ensembl_write_assemblies(connections, query, example_file, example_query, output,
                             db_url=None):
    """
    This will query the ensembl databases in the connections file and write the
    output to the given file.
    """
    assemblies.write(connections, query, example_file, example_query, db_url, output)


@cli.command('coordinate-systems')
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql', type=click.File('r'))
@click.argument('output', default='coordinate_systems.csv', type=click.File('w'))
def ensembl_coordinates(connections, query, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    coordinate_systems.write(connections, query, output)


@cli.command('karyotypes')
@click.argument('output', default='karyotypes.csv', type=click.File('w'))
@click.argument('species', nargs=-1)
def ensembl_write_karyotypes(output, species):
    """
    Fetch all the karyotype information from all Ensembl species. This will use
    the Ensembl API to fetch the data and write to the given output file.
    """
    if not species:
        species = None
    else:
        species = set(species)
    karyotypes.write(output, species=species)


@cli.command('proteins')
@click.argument('connections', default='databases.json', type=click.File('r'))
@click.argument('query', default='query.sql-', type=click.File('r'))
@click.argument('output', default='proteins.csv', type=click.File('w'))
def ensembl_proteins_cmd(connections, query, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    proteins.write(connections, query, output)


@cli.command('compara')
@click.argument('fasta', default='-', type=click.File('rb'))
@click.argument('output', default='compara.csv', type=click.File('wb'))
def ensembl_compara(fasta, output):
    """
    Parse the FASTA file of Ensembl compara data. This will produce a CSV file
    for import into the database that tracks what ensembl transcripts are
    homologous.
    """
    compara.write(fasta, output)
