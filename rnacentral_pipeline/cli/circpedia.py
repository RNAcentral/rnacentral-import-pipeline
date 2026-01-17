# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.circpedia import parser
from rnacentral_pipeline.writers import entry_writer


@click.group("circpedia")
def cli():
    """
    Commands for parsing CIRCpedia V3 circular RNA data.
    """


@cli.command("parse")
@click.argument("taxonomy", type=click.Path(exists=True))
@click.argument("csv_file", type=click.Path(exists=True))
@click.argument(
    "output",
    default=".",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
)
@click.option(
    "--assembly",
    default=None,
    help="Genome assembly ID (e.g., GRCh38, GRCm39)",
)
def parse_circpedia(taxonomy, csv_file, output, assembly):
    """
    Parse CIRCpedia V3 CSV data file.

    This command parses circular RNA data from CIRCpedia V3 and converts it
    to RNAcentral import format.

    Arguments:

        TAXONOMY: Path to taxonomy database file (context.db)

        CSV_FILE: Path to CIRCpedia CSV data file

        OUTPUT: Output directory for CSV files (default: current directory)

    The CSV file should contain columns:
        - circid: CIRCpedia circular RNA ID
        - species: Species name
        - location: Genomic location (chr:start-end)
        - strand: Strand (+ or -)
        - gene: Host gene name (optional)
        - sequence: RNA sequence (optional)
        - exon_positions: Exon coordinates (optional)
        - fpm: Expression level (optional)
        - cell_line: Cell line/tissue (optional)
    """
    entries = parser.parse(csv_file, Path(taxonomy), assembly_id=assembly)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
