# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from pathlib import Path

import click

from rnacentral_pipeline.writers import entry_writer

from . import parser


@click.group("tarbase")
def cli():
    """
    Commands for parsing Tarbase V9 data.
    """


@cli.command("parse")
@click.argument("tsv_file", type=click.File("r"))
@click.argument("fasta", type=click.Path())
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
def process_tsv(json_file, fasta, output):
    """
    Parse the TSV file TarBase provides and a fasta file of miRNA sequences to
    produce importable CSV files. This requires the files be sorted by the
    mirna_name column.
    """
    entries = parser.parse(json_file, Path(fasta))
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
