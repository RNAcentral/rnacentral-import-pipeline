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

import json

import click


@click.group('pdb')
def cli():
    """
    Commands for dealing with the fetching PDB data.
    """
    pass


@cli.command('data')
@click.argument('output', default='pdb.json', type=click.File('wb'))
def pdb_group_data(output):
    pdb_ids = pdb.rna_containing_pdb_ids()
    json.dump(pdb.ncrna_chain_descriptions(pdb_ids), output)


@cli.command('extra')
@click.argument('output', default='pdb-extra.json', type=click.File('wb'))
def pdb_group_extra(output):
    pdb_ids = pdb.rna_containing_pdb_ids()
    json.dump(pdb.reference_mapping(pdb_ids), output)
