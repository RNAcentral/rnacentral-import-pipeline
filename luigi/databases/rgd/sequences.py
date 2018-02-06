# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import gzip
import shutil
import tempfile

from Bio import SeqIO
from pybedtools import BedTool


def merge_gff(filenames, handle):
    header, rest = filenames
    with gzip.open(header, 'rb') as raw:
        shutil.copyfileobj(raw, handle, 1024*1024*10)

    for filename in rest:
        with gzip.open(filename, 'rb') as raw:
            for line in raw:
                if line.startswith('#'):
                    continue
                handle.write(line)
    handle.seek(0)
    return BedTool(handle)


def merge_fasta_files(filenames, handle):
    for filename in filenames:
        with gzip.open(filename, 'rb') as raw:
            shutil.copyfileobj(raw, handle, 1024*1024*10)


def extract(gff_files, fasta_files, filename):
    with tempfile.NamedTemporaryFile() as tmp_gff, \
            tempfile.NamedTemporaryFile() as fasta:
        gff = merge_gff(gff_files, tmp_gff)
        merge_fasta_files(fasta_files, fasta)
        result = gff.sequences(fi=fasta, s=True, fo=filename)
        return SeqIO.index(result.seqfn, 'fasta')
