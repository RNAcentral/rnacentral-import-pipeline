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

from rnacentral_pipeline.databases.helpers import publications


@click.group("europepmc")
def cli():
    """
    This group of commands deal with handling publication data.
    """
    pass


@cli.command('index-xml')
@click.argument('directory', default='out', type=click.Path())
@click.argument('db', default='references.db', type=click.Path())
def index_xml(directory, db):
    """
    Index a list of XML files of publication information.
    """
    publications.index_xml_directory(directory, db)


@cli.command('lookup')
@click.option('--column', default=0)
@click.option('--allow-fallback/--no-allow-feedback', default=False)
@click.argument('db', default='references.db', type=click.Path())
@click.argument('ids', default='ref_ids.csv', type=click.File('r'))
@click.argument(
    'output', 
    default='references.csv', 
    type=click.File('w'))
def lookup(db, ids, output, column=0, allow_fallback=False):
    """
    Use the database index file to lookup all reference information for all xml
    files in the given directory.
    """
    publications.write_file_lookup(db, ids, output, column=column, allow_fallback=allow_fallback)
