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

import tempfile
from ftplib import FTP
import subprocess as sp
from contextlib import contextmanager

from Bio import SeqIO
from Bio import Entrez

from boltons import iterutils

from rnacentral_pipeline.utils import pickle_stream

from . import helpers

BATCH_SIZE = 200


@contextmanager
def raw():
    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()
    ftp.cwd('gene/DATA')
    with tempfile.NamedTemporaryFile(suffix='.gz') as gz:
        ftp.retrbinary('RETR gene_info.gz', gz.write)
        gz.flush()
        ftp.quit()
        with tempfile.NamedTemporaryFile(mode='w+') as tmp:
            sp.run(['gzip', '-cd', gz.name], check=True, stdout=tmp)
            tmp.flush()
            tmp.seek(0)
            yield tmp


def ncrnas():
    with raw() as handle:
        for ncrna in helpers.ncrnas(handle):
            yield ncrna


def sequences(ncrna_ids):
    return {}


def data():
    found = 0
    total = 0
    batches = iterutils.chunked_iter(ncrnas(), BATCH_SIZE)
    for batch in batches:
        seqs = sequences(batch)
        for ncrna in batch:
            total += 1
            gene_id = helpers.gene_id(ncrna)
            if gene_id not in seqs:
                LOGGER.warn("No sequence found for %s" % gene_id)
                continue
            found += 1
            ncrna['sequence'] = seqs[gene_id]
            yield ncrna

    if float(found) / float(total) <= 0.95:
        raise ValueError("Not enough sequences found")


def write(output):
    return pickle_stream(data(), output)
