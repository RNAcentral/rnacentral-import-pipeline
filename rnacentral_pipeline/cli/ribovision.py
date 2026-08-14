# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from pathlib import Path

import click
from Bio import SeqIO

from rnacentral_pipeline.databases.ribovision import helpers, parser
from rnacentral_pipeline.writers import entry_writer


@click.group("ribovision")
def cli():
    """
    Commands for parsing ribovision data
    """
    pass


@cli.command("parse")
@click.argument("metadata_file", type=click.File("r"))
@click.argument("sequence_file", type=click.Path(dir_okay=False))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_ribovision(metadata_file, sequence_file, output):
    """
    Build entries from the r2dt_models metadata and the sequences R2DT ships.
    """
    entries = parser.parse(metadata_file, Path(sequence_file))
    with entry_writer(Path(output)) as writer:
        writer.write(entries)


@cli.command("r2dt-to-fasta")
@click.argument("directory", type=click.Path())
@click.argument("output", type=click.File("w"))
def generate_r2dt_fasta(directory, output):
    """
    Collect the LSU and SSU model sequences from an R2DT data directory.
    """
    entries = helpers.fasta_entries(Path(directory))
    SeqIO.write(entries, output, "fasta")
