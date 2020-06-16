#!/usr/bin/env python3

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

import re
import json

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

EASEL_PATTERN = re.compile(r'^[ACGUN]$')


def as_record(entry):
    description = entry.get('description', '') or ''

    return SeqRecord(
        Seq(entry['sequence']),
        id=entry['id'],
        description=description,
    )


def select_easel(sequence):
    return re.match(EASEL_PATTERN, entry['sequence'])


def sequences(handle, only_valid_easel=False):
    """
    Parse each line and create a generator of SeqRecords to write.
    """

    data = map(json.loads, handle)
    if only_valid_easel:
        data = filter(select_easel, data)
    return map(as_record, data)


@click.command()
@click.option('--only-valid-easel', default=False)
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('w'))
def cli(tsv_file, output=None, only_valid_easel=False):
    """
    Turn a file with JSON lines into a FASTA file. Each line must be an object
    with an id and sequence entry and an optional description entry.
    """
    seqs = sequences(tsv_file, only_valid_easel=only_valid_easel)
    SeqIO.write(seqs, output, "fasta")


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
