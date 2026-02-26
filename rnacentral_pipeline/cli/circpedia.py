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
@click.argument("annotation_file", type=click.Path(exists=True))
@click.argument("fasta_file", type=click.Path(exists=True))
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
def parse_circpedia(annotation_file, fasta_file, output, assembly):
    """
    Parse CIRCpedia V3 annotation and sequence files.

    This command parses circular RNA data from CIRCpedia V3 and converts it
    to RNAcentral import format.

    Arguments:

        ANNOTATION_FILE: Path to CIRCpedia TSV annotation file

        FASTA_FILE: Path to CIRCpedia FASTA sequence file

        OUTPUT: Output directory for CSV files (default: current directory)

    The annotation file should be tab-delimited with columns:
        - circID: CIRCpedia circular RNA ID
        - species: Species name
        - Location: Genomic location with strand (e.g., V:15874634-15876408(-))
        - gene_Ensembl: Ensembl gene ID (optional)
        - gene_Refseq: RefSeq gene ID (optional)
        - circname: circRNA name (optional)
        - subcell_location: Subcellular localization (optional)
        - editing_site: Editing sites (optional)
        - DIS3_signal: DIS3 degradation signals (optional)
        - Orthology: Orthology information (optional)
        - TGS: Third-generation sequencing support (optional)

    The FASTA file should contain sequences with headers matching circIDs.
    """
    entries = parser.parse(annotation_file, fasta_file, assembly_id=assembly)
    with entry_writer(Path(output)) as writer:
        writer.write(entries)
