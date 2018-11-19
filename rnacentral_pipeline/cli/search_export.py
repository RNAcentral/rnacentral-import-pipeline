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

from rnacentral_pipeline.rnacentral.search_export import exporter as search


@click.group('search-export')
def cli():
    """
    A group of commands dealing with our search export pipeline.
    """
    pass


@cli.command('as-xml')
@click.argument('raw_file', type=click.File('rb'))
@click.argument('metadata_file', type=click.File('rb'))
@click.argument('xml_file', type=click.File('wb'))
@click.argument('count_file', type=click.File('wb'), default='count')
def search_export_xml(raw_file, metadata_file, xml_file, count_file=None):
    """
    This will parse a file with one JSON object per line to produce XML
    formatted data that is used as input to the search team. Additionally, this
    produces a count file which contains the number of entries in the XML file.
    This is needed for building the release_note.txt file.
    """
    search.as_xml(raw_file, metadata_file, xml_file, count_file)


@cli.command('release-note')
@click.argument('release', type=str)
@click.argument('output', type=click.File('wb'))
@click.argument('count_files', nargs=-1, type=click.File('rb'))
def search_export_note(release, output, count_files):
    """
    This will create the release_note.txt file that is needed for the search
    export.
    """
    search.release_note(output, release, count_files)
