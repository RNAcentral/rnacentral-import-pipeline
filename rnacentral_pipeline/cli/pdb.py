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

import logging

import click

from rnacentral_pipeline.databases.pdb import fetch
from rnacentral_pipeline.databases.pdb import parser
from rnacentral_pipeline.writers import write_entries

LOGGER = logging.getLogger(__name__)


@click.group('pdb')
def cli():
    """
    Commands for dealing with the fetching PDB data.
    """
    pass


@cli.command("generate")
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
    chain_info = fetch.rna_chains()
    references = {}
    try:
        references = fetch.references(chain_info)
    except Exception:
        LOGGER.info("Failed to get extra references")
    write_entries(parser.parse, output, chain_info, references)
