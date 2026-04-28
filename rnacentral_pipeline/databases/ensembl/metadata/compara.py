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

import csv
import hashlib
from pathlib import Path

from Bio import AlignIO
from six.moves import cStringIO as StringIO

from rnacentral_pipeline import schemas
from rnacentral_pipeline.output_format import is_parquet
from rnacentral_pipeline.parquet_writers import parquet_writer


def alignments(fasta):
    """
    Parse the given fasta file and produce an iterable of aligned sequences for
    each chunk of aligned sequences.
    """

    buff = StringIO()
    for line in fasta:
        line = str(line)
        if line.startswith("//"):
            buff.seek(0)
            yield AlignIO.read(buff, "fasta")
            buff = StringIO()
        else:
            buff.write(line)


def data(fasta):
    """
    Produce an iterable for each chunk of aligned sequences. This will add a
    uuid to each chunk to track which sequences are homologous.
    """

    for alignment in alignments(fasta):
        ids = []
        for entry in alignment:
            ids.append(entry.id)
        ids.sort()

        unique = hashlib.sha256()
        for eid in ids:
            unique.update(eid.encode("utf-8"))
        yield (unique.hexdigest(), ids)


def _rows(fasta):
    for key, ids in data(fasta):
        for eid in ids:
            yield (key, eid)


def write(fasta, output):
    """
    Write the homology data in the given fasta file to the given output file.
    """

    if isinstance(output, (str, Path)):
        path = Path(output)
        if is_parquet():
            with parquet_writer(path, schemas.COMPARA) as writer:
                for row in _rows(fasta):
                    writer.writerow(row)
            return
        with path.open("w") as handle:
            writer = csv.writer(handle)
            for row in _rows(fasta):
                writer.writerow(list(row))
        return

    writer = csv.writer(output)
    for row in _rows(fasta):
        writer.writerow(list(row))
