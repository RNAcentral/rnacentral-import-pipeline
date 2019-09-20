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

import io
from ftplib import FTP

from Bio import SeqIO
from Bio import Entrez

from boltons import iterutils

from rnacentral_pipeline.utils import pickle_stream

from . import helpers

BATCH_SIZE = 200


def ncrnas():
    ftp = FTP('ftp://ftp.ncbi.nih.gov/gene/')
    ftp.login()
    handle = io.StringIO()
    ftp.retrbinary('RETR %s' % filename, handle.write)
    ftp.quit()
    handle.seek(0)
    return helpers.ncrnas(handle)


def sequences(ncrna_ids):
    return {}
    # Entrez.esearch('


def fetch():
    found = {helpers.gene_id(ncrna): ncrna for ncrna in ncrnas()}
    batches = iterutils.chunked_iter(found.values(), BATCH_SIZE)
    for batch in batches:
        seqs = sequences(batch)
        for gene_id, ncrna in batch:
            if gene_id not in seqs:
                LOGGER.warn("No sequence found for %s" % gene_id)
                continue
            ncrna['sequence'] = seqs[gene_id]
            yield ncrna


def write(output):
    return pickle_stream(fetch(), handle)
