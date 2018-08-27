#!/usr/bin/env python

import json

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def sequences(handle):
    for line in handle:
        data = json.loads(line)
        description = data.get('description', None)
        if description:
            description = unicode(description)

        yield SeqRecord(
            Seq(data['sequence']),
            id=data['id'].encode('utf-8'),
            description=description,
        )


@click.command()
@click.argument('tsv_file', type=click.File('rb'))
@click.argument('output', default='-', type=click.File('wb'))
def cli(tsv_file, output=None):
    SeqIO.write(sequences(tsv_file), output, "fasta")


if __name__ == '__main__':
    cli()  # pylint: disable=no-value-for-parameter
