#!/usr/bin/env python

import csv

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sequences(handle):
    for row in csv.reader(tsv_file, sep='\t'):
        yield SeqRecord(Seq(row[1]), id=row[0])


@click.command()
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', type=click.File('wb'))
def cli(tsv_file, output):
    SeqIO.write(sequences(tsv_file), output, "fasta")
