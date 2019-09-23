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

BATCH_SIZE = 300

Entrez.email = 'rnacentral@gmail.com'


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


def find_nucleotide_id(entry):
    seq_id = None
    cur_id = entry['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
    for locus in entry['Entrezgene_locus']:
        for product in locus.get('Gene-commentary_products', []):
            if 'Gene-commentary_rna' in product \
                    and 'Gene-commentary_accession' in product:
                if seq_id is not None:
                    raise ValueError("Found duplicate sequence id for: %s" % cur_id)
                seq_id = product['Gene-commentary_accession']
    return seq_id


def find_genome_id(entry):
    seq_id = None
    start = None
    stop = None
    cur_id = entry['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
    for locus in entry['Entrezgene_locus']:
        pass
    if seq_id is None or start is None or stop is None:
        return None
    return (seq_id, start, stop)


def lookup_by_nt(mapping):
    handle = Entrez.efetch(
        db='nucleotide', 
        id=','.join(mapping.keys()), 
        rettype="gb",
        retmode="text",
    )
    data = {}
    for sequence in SeqIO.parse(handle, 'genbank'):
        # The ID has a version ID, which we do not want, while name does not.
        gene_id = mapping[sequence.name]
        data[gene_id] = sequence
    handle.close()
    assert len(data) == len(mapping)
    return data


def lookup_by_genome(mapping):
    return {}


def sequences(batch):
    ids = [ncrna['GeneID'] for ncrna in batch]
    handle = Entrez.efetch(db='gene', id=ids, retmode='xml')
    ncrna_ids = {}
    genomic_ids = {}
    for entry in Entrez.read(handle):
        cur_id = entry['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']
        nt_id = find_nucleotide_id(entry)
        if not nt_id:
            genome_id = find_genome_id(entry)

        genome_id = None
        if nt_id:
            ncrna_ids[nt_id] = cur_id
        elif genome_id:
            genomic_ids[genome_id] = cur_id
        else:
            print("Could not find sequence id for: %s" % cur_id)
            # raise ValueError("Could not find sequence id for: %s" % cur_id)

    data = {}
    data.update(lookup_by_nt(ncrna_ids))
    data.update(lookup_by_genome(genomic_ids))
    return data


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
