#!/usr/bin/env python

import csv

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sequences(handle):
    for line in handle:
        row = line.split('\t')
        description = ''
        if len(row) > 2:
            description = row[2]
        yield SeqRecord(Seq(row[1]), id=row[0], description=description)


@click.command()
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def cli(tsv_file, output=None):
    SeqIO.write(sequences(tsv_file), output, "fasta")


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
