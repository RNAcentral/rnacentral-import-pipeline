#!/usr/bin/env python

import csv

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sequences(handle):
    for row in csv.reader(handle, delimiter='\t'):
        yield SeqRecord(Seq(row[1]), id=row[0], description='')


@click.command()
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def cli(tsv_file, output=None):
    SeqIO.write(sequences(tsv_file), output, "fasta")


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
