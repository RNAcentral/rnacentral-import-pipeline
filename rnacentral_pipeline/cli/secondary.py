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

from rnacentral_pipeline.rnacentral import secondary


@click.group('secondary')
def cli():
    """
    A group of commands for parsing data from secondary structures into an
    importable format.
    """
    pass


@cli.command('process-svgs')
@click.argument('directories', nargs=-1, type=click.Path(
    writable=True,
    dir_okay=True,
    file_okay=False,
))
@click.argument('output', type=click.File('w'))
def process_svgs(directories, output):
    """
    Process all SVG secondary structures in the given directory and produce a
    single data file that can be imported into the database.
    """
    secondary.write_all(directories, output)