# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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
import pickle

import click

from rnacentral_pipeline.databases import pdb
from rnacentral_pipeline.databases.pdb import parser
from rnacentral_pipeline.writers import write_entries


@click.group('pdb')
def cli():
    """
    Commands for dealing with the fetching PDB data.
    """
    pass


@cli.group('fetch')
def fetch():
    """
    Commands for fetch PDB(e) data.
    """


@fetch.command('data')
@click.argument('output', default='pdb.json', type=click.File('w'))
@click.argument('pdb_ids', nargs=-1)
def pdb_group_data(output, pdb_ids=None):
    json.dump(pdb.all_rna_chains(), output)


@fetch.command('extra')
@click.argument('chains', type=click.File('r'))
@click.argument('output', default='pdb-extra.json', type=click.File('wb'))
@click.argument('pdb_ids', nargs=-1)
def pdb_group_extra(chains, output, pdb_ids=None):
    data = json.load(chains)
    pickle.dump(pdb.references(data), output)


@cli.command("parse")
@click.argument("pdb_data", default="pdb.json", type=click.File("rb"))
@click.argument("extra", default="pdb-extra.json", type=click.File("rb"))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False,),
)
def process_pdb(pdb_data, extra, output):
    """
    This will fetch and parse all sequence data from PDBe to produce the csv
    files we import.
    """
    write_entries(parser.parse, output, pdb_data, extra)
