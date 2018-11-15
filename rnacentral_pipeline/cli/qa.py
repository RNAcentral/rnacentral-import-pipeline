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

from rnacentral_pipeline.databases import rfam
from rnacentral_pipeline.databases import dfam


@click.group("qa")
def cli():
    """
    This group of commands deal with QA work.
    """
    pass


@cli.command('rfam')
@click.argument('tblout', default='-', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def process_tblout(tblout, output):
    """
    Process a table out file and create a CSV for importing into our database.
    This will overwrite the given file.
    """
    rfam.infernal_results.as_csv(tblout, output)


@cli.command('dfam')
@click.argument('data', default='-', type=click.File('r'))
@click.argument('output', default='-', type=click.File('w'))
def process_dfam(data, output):
    """
    This will parse the results from a dfamscan.pl call and turn the results
    into a CSV file. This assumes the results are in the space delimited
    format.
    """
    dfam.as_csv(data, output)
