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

from rnacentral_pipeline.databases.europepmc import xml


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
    xml.index_directory(directory, db)


# @cli.command('query')
# @click.option('--column', default=0)
# @click.option('--ignore-missing/--no-ignore-missing', default=True)
# @click.argument('references', default='-', type=click.File('r'))
# @click.argument('output', default='-', type=click.File('w'))
# def query(references, output):
#     """
#     Query EuropePMC API for all references in the references file. This will be
#     subject to rate limits to ensure we do not overload the API so it will take
#     a long time for large files. 
#     """
#     fetch.write_file_lookup(
#         references, 
#         output,
#         column=column, 
#         ignore_missing=ignore_missing,
#     )


@cli.command('lookup')
@click.option('--column', default=0)
@click.option('--allow-fallback/--no-allow-fallback', default=False)
@click.option('--ignore-missing/--no-ignore-missing', default=True)
@click.argument('db', default='references.db', type=click.Path())
@click.argument('ids', default='ref_ids.csv', type=click.File('r'))
@click.argument(
    'output', 
    default='references.csv', 
    type=click.File('w'))
def lookup(db, ids, output, column=0, allow_fallback=False, ignore_missing=True):
    """
    Use the database index file to lookup all reference information for all xml
    files in the given directory.
    """
    xml.write_file_lookup(
        db, 
        ids, 
        output, 
        column=column, 
        allow_fallback=allow_fallback,
        ignore_missing=ignore_missing,
    )
