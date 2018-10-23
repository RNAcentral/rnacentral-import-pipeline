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

from rnacentral_pipeline.databases.ensembl import tasks as ensembl_tasks
from rnacentral_pipeline.databases.ensembl import proteins as ensembl_proteins
from rnacentral_pipeline.databases.ensembl import assemblies as ensembl_assemblies
from rnacentral_pipeline.databases.ensembl import coordinate_systems as ensembl_coords


@click.group('ensembl')
def cli():
    """
    This is a set of commands for dealing with processing protein information.
    We don't have much in the way of protein summary but sometimes we do need a
    little for display.
    """
    pass


@cli.command('proteins')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_proteins_cmd(filename, output):
    """
    This will process the ensembl protein information files. This assumes the
    file is sorted.
    """
    ensembl_proteins.from_file(filename, output)


@cli.command('coordinate-systems')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_coordinates(filename, output):
    """
    Turn the tsv from the ensembl query into a csv that can be imported into
    the database.
    """
    ensembl_coords.from_file(filename, output)


@cli.command('select-tasks')
@click.argument('filename', type=click.File('rb'))
@click.argument('possible', type=click.File('rb'))
@click.argument('done', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_select_tasks(filename, possible, done, output):
    """
    This will select which database names to import. It basically checsk to see
    which is the newest assuming Ensembl's standard naming schema and see which
    we have not already imported.
    """
    ensembl_tasks.write(filename, possible, done, output)


@cli.command('assemblies')
@click.argument('filename', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def ensembl_write_assemblies(filename, output):
    ensembl_assemblies.write(filename, output)
