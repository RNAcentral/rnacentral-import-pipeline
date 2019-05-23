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

from rnacentral_pipeline.databases.ols import helpers as ols
from rnacentral_pipeline.writers import write_ontology_annotations as onto_writer
from rnacentral_pipeline.databases.quickgo import parser as quickgo


@click.group('ontologies')
def cli():
    """
    Entry point for fetching ontology data.
    """
    pass


@cli.command('quickgo')
@click.argument('raw_data', type=click.File('r'))
@click.argument('output', default='.', type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
def ontologies_quickgo(raw_data, output):
    """
    This will process a quickgo file and output files into the given directory.
    """
    onto_writer(quickgo.parser, output, raw_data)


@cli.command('lookup-terms')
@click.argument('terms', type=click.File('r'))
@click.argument('output', type=click.File('w'))
def ontologies_lookup_terms(terms, output):
    ols.lookup_terms(terms, output)
