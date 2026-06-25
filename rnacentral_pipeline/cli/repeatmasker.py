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

from rnacentral_pipeline.repeatmasker import parser
from rnacentral_pipeline.repeatmasker.data import RepeatmaskerWriter
from rnacentral_pipeline.writers import build


@click.group("repeatmasker")
def cli():
    """
    Commands for parsing RepeatMasker output.
    """
    pass


@cli.command("parse")
@click.argument("sequences", type=click.Path())
@click.argument("repeatmasker_out", type=click.Path())
@click.argument("output", type=click.Path())
def parse(sequences, repeatmasker_out, output):
    """
    Parse a RepeatMasker .out file against the scanned SEQUENCES FASTA, writing
    results.csv (per-URS presence) and features.csv (per-region coordinates).
    """
    data = parser.parse(Path(sequences), Path(repeatmasker_out))
    with build(RepeatmaskerWriter, Path(output)) as wtr:
        wtr.write(data)
