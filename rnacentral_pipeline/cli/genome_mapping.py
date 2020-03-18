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
from rnacentral_pipeline.rnacentral.genome_mapping import attempted


@click.group('genome-mapping')
def cli():
    """
    This group of commands deals with figuring out what data to map as well as
    parsing the result into a format for loading.
    """
    pass


@cli.group('blat')
def hits():
    """
    A series of commands for working with blat hits.
    """
    pass


@hits.command('serialize')
@click.argument('assembly_id')
@click.argument('hits', default='-', type=click.File('r'))
@click.argument('output', default='-', type=click.File('wb', lazy=False))
def hits_json(assembly_id, hits, output):
    """
    Serialize the PSL file into something that python can later process. This is
    a lossy operation but keeps everything needed for selecting later. This
    exists so we can do mulitple select steps and still merge the results.
    """
    blat.as_pickle(assembly_id, hits, output)


@hits.command('as-importable')
@click.argument('hits', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('w', lazy=False))
def as_importable(hits, output):
    """
    Convert a json-line file into a CSV that can be used for import by pgloader.
    This is lossy as it only keeps the things needed for the database.
    """
    blat.write_importable(hits, output)


@hits.command('select')
@click.option('--sort', is_flag=True, default=False)
@click.argument('hits', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb', lazy=False))
def select_hits(hits, output, sort=False):
    """
    Parse a JSON-line file and select the best hits in the file. The best hits
    are written to the output file. This assumes the file is sorted by
    urs_taxid unless --sort is given in which case the data is sorted in memory.
    That may be very expensive.
    """
    blat.select_pickle(hits, output, sort=sort)


@cli.command('url-for')
@click.option('--host', default='ensembl')
@click.argument('species')
@click.argument('assembly_id')
@click.argument('output', default='-', type=click.File('w'))
def find_remote_url(species, assembly_id, output, host=None):
    """
    Determine the remote URL to fetch a the genome for a given species/assembly.
    The url is written to the output file and may include '*'.
    """
    url = urls.url_for(species, assembly_id, host=host)
    output.write(url)


@cli.command('urls-for')
@click.argument('filename', default='-', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def find_remote_urls(filename, output):
    """
    Determine the remote URL to fetch a the genomes for all entries in a file,
    where the file is a csv of species,assembly. The urls is written to the
    output file and may include '*'.
    """
    urls.write_urls_for(filename, output)


@cli.command('create-attempted')
@click.argument('filename', type=click.File('r'))
@click.argument('taxid')
@click.argument('assembly_id')
@click.argument('output', type=click.File('w'))
def parse_attempted_sequences(filename, taxid, assembly_id, output):
    attempted.write(filename, taxid, assembly_id, output)
