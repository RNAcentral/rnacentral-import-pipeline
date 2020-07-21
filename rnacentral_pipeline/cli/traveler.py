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

from rnacentral_pipeline.rnacentral import traveler
from rnacentral_pipeline.rnacentral import attempted

from rnacentral_pipeline.databases.crw.traveler import models as crw
from rnacentral_pipeline.databases.ribovision.traveler import models as ribovision
from rnacentral_pipeline.databases.gtrnadb.traveler import models as gtrnadb

from . import parameters as params


@click.group('traveler')
def cli():
    """
    A group of commands for parsing data from secondary structures into an
    importable format.
    """
    pass


@cli.command('process-svgs')
@click.option('--allow-missing', is_flag=True, default=False)
@click.argument('model_info', type=click.File('r'))
@click.argument('directory', type=click.Path(
    writable=False,
    dir_okay=True,
    file_okay=False,
))
@click.argument('output', type=click.File('w'))
def process_svgs(model_info, directory, output, allow_missing=False):
    """
    Process all SVG secondary structures in the given directory and produce a
    single data file that can be imported into the database.
    """
    traveler.write(model_info, directory, output, allow_missing=allow_missing)


@cli.group('model-info')
def model_info():
    """
    Commands for parsing and generating data files we can import into the
    database as model info files.
    """
    pass


@model_info.command('crw')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def crw_model_info(filename, output):
    """
    Parse the CRW metadata file and produce 
    """
    crw.write(filename, output)


@model_info.command('ribovision')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def ribovision_model_info(filename, output):
    """
    Parse the metadata.tsv file from auto-traveler for Ribovision models to
    produce something we can put in our database.
    """
    ribovision.write(filename, output)


@model_info.command('gtrnadb')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def gtrnadb_model_info(filename, output):
    """
    Parse the metadata.tsv file from auto-traveler for gtrnadb models to
    produce something we can put in our database.
    """
    gtrnadb.write(filename, output)


@cli.command('create-attempted')
@click.argument('filename', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def traveler_create_attempted(filename, output):
    attempted.traveler(filename, output)


@cli.command('publish')
@click.option('--allow-missing', is_flag=True, default=False)
@click.argument('model_info', type=click.File('r'))
@click.argument('directory', type=click.Path(
    writable=False,
    dir_okay=True,
    file_okay=False,
))
@click.argument('output', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def traveler_publish(model_info, directory, output, allow_missing):
    traveler.publish(model_info, directory, output, allow_missing=allow_missing)
