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

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import publications as pubs

from .data import KnownMapper

RNA_TYPE_MAPPING = {
    'RNA, Y': 'Y_RNA',
    'RNA, cluster': 'other',
    'RNA, long non-coding': 'lncRNA',
    'RNA, micro': 'precursor_RNA',
    'RNA, misc': 'misc_RNA',
    'RNA, ribosomal': 'rRNA',
    'RNA, small cytoplasmic': 'scRNA',
    'RNA, small nuclear': 'snRNA',
    'RNA, small nucleolar': 'snoRNA',
    'RNA, transfer': 'tRNA',
    'RNA, vault': 'vault_RNA',
}


def known_refseq(dbconf):
    pass


def known_gtrnadb(dbconf):
    pass


def known_ensembl(dbconf):
    pass


def known(dbconf):
    mapper = KnownMapper()
    mapper.store_all(known_refseq(dbconf))
    mapper.store_all(known_gtrnadb(dbconf))
    mapper.store_all(known_ensembl(dbconf))
    return mapper


def gene_symbol(raw):
    if raw['prev_symbol']:
        return ';'.join(raw['prev_symbol'])
    return None


def organelle(raw):
    if raw['location'] == 'mitochondria':
        return 'Mitochondrion'
    return None


def description(raw):
    gene_description = ''
    if raw['symbol']:
        gene_description = ' (%s)' % raw['symbol']
    return 'Homo sapiens (human)' + raw['name'] + gene_description


def rna_type(raw):
    return RNA_TYPE_MAPPING[raw['locus_type']]


def references(raw):
    return [pubs.reference(pmid) for pmid in raw.get('pubmed_id', [])]


def as_entry(raw):
    return data.Entry(
        id=raw['hgnc_id'],
        rna_type=rna_type(raw),

        database='HGNC',
        external_id=raw['symbol'],
        species='Homo sapiens',
        organelle=organelle(raw),
        gene=raw['symbol'],
        gene_synonym=gene_symbol(raw),
        description=description(raw),
        references=references(raw),
    )
