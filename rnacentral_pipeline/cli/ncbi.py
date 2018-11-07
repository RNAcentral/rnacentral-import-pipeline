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

from rnacentral_pipeline.databases.ncbi import taxonomy


@click.group("ncbi")
def cli():
    """
    This group of commands deal with processing NCBI data.
    """
    pass


@cli.command("taxonomy")
@click.argument('lineage', default='fulllineage.dmp', type=click.File('r'))
@click.argument('names', default='names.dmp', type=click.File('r'))
@click.argument('merged', default='merged.dmp', type=click.File('r'))
@click.argument('output', default='taxonomy.csv', type=click.File('w'))
def parse_taxonomy(lineage, names, merged, output):
    taxonomy.write(lineage, names, merged, output)
