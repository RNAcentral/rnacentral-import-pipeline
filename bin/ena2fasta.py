#!/usr/bin/env python

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

from os.path import getsize

import click

from Bio import SeqIO


def records(filename):
    for record in SeqIO.parse(filename, "embl"):
        if 'metagenome' in lineage(record).lower():
            yield record


@click.command()
@click.argument('filename', type=click.File('r'))
@click.argument('output', type=click.File('w'))
def main(filename, output):
    """
    Convert a ENA EMBL file into a fasta file suitable for ribotyper analysis.
    """
    SeqIO.write(sequences(filename), output, "fasta")
    output.flush()

    if not getsize(output.name):
        os.remove(output.name)


if __name__ == '__main__':
    main()
