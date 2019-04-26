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

from rnacentral_pipeline.rnacentral.genome_mapping import urls
from rnacentral_pipeline.rnacentral.genome_mapping import blat


@click.group('genome-mapping')
def cli():
    """
    This group of commands deals with figuring out what data to map as well as
    parsing the result into a format for loading.
    """
    pass


@cli.command('select-hits')
@click.argument('assembly_id')
@click.argument('hits', default='-', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def select_hits(assembly_id, hits, output):
    blat.write_selected(assembly_id, hits, output)


@cli.command('url-for')
@click.option('--host', default='ensembl')
@click.argument('species')
@click.argument('assembly_id')
@click.argument('output', default='-', type=click.File('w'))
def find_remote_url(species, assembly_id, output, host=None):
    url = urls.url_for(species, assembly_id, host=host)
    output.write(url)


@cli.command('urls-for')
@click.argument('filename', default='-', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def find_remote_urls(filename, output):
    urls.write_urls_for(filename, output)
