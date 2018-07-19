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

import re

from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs

from .data import KnownMapper
from .data import PartialEntry


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

ONE_TO_THREE = {
    'A': 'Ala',
    'C': 'Cys',
    'D': 'Asp',
    'E': 'Glu',
    'F': 'Phe',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'K': 'Lys',
    'L': 'Leu',
    'M': 'Met',
    'N': 'Asn',
    'P': 'Pro',
    'Q': 'Gln',
    'R': 'Arg',
    'S': 'Ser',
    'T': 'Thr',
    'U': 'SeC',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr',
    'X': 'iMet',
    'SUP': 'Sup',
}


def known_refseq(dbconf):
    pass


def known_gtrnadb(dbconf):
    pass


def known_ensembl(dbconf):
    pass


def known_ensembl_sequences(dbconf):
    pass


def known(dbconf):
    mapper = KnownMapper()
    mapper.store_all(known_refseq(dbconf))
    mapper.store_all(known_gtrnadb(dbconf))
    mapper.store_all(known_ensembl(dbconf))
    return mapper


def gene_symbols(raw):
    """
    Get a list of previous symbols this gene has been known as.
    """

    if raw['prev_symbol']:
        return raw['prev_symbol']
    return []


def organelle(raw):
    """
    Determine what organelle the entry is in. This can detect if a sequence
    comes from a mitochondria.
    """

    if raw['location'] == 'mitochondria':
        return 'Mitochondrion'
    return None


def description(raw):
    """
    Build a description for RNAcentral out of the given name and gene symbol.
    """

    gene_description = ''
    if raw['symbol']:
        gene_description = ' (%s)' % raw['symbol']
    return 'Homo sapiens (human)' + raw['name'] + gene_description


def rna_type(raw):
    """
    Determine the standardized RNA types from the given locus type.
    """
    return RNA_TYPE_MAPPING[raw['locus_type']]


def references(raw):
    """
    Compute all data needed for the references that are present in the given
    HGNC reference.
    """
    return [pubs.reference(pmid) for pmid in raw.get('pubmed_id', [])]


def taxid(_):
    """
    Returns 9606 every time.
    """
    return 9606


def lineage(raw):
    """
    Compute the lineage, will use an API the first time.
    """
    return phy.lineage(taxid(raw))


def common_name(raw):
    """
    Compute a common name, will use an API the first time.
    """
    return phy.common_name(taxid(raw))


def species(raw):
    """
    Compute the species, will use an API the first time.
    """
    return phy.species(taxid(raw))


def gtrnadb_id(raw):

    prefix = None
    accession = raw['symbol']
    match = re.match(r'TR(\S+)-(\S{3})(\d+-\d+)', accession)
    if match:
        prefix = 'tRNA'

    # nuclear-encoded mitochondrial tRNAs
    match = re.match(r'NMTR(\S+)-(\S{3})(\d+-\d+)', accession)
    if match:
        prefix = 'nmt-tRNA'

    if not prefix:
        return None

    return '{prefix}-{anticodon}-{word}-{range}'.format(
        prefix=prefix,
        anticodon=ONE_TO_THREE[match.group(1)],
        word=match.group(2),
        range=match.group(3),
    )


def ensembl_id(raw):
    pass


def refseq_id(raw):
    return raw.get('refseq_accession', None)


def xref_data(raw):
    """
    Compute xrefs from HGNC. This will only contain the xrefs that we can map
    and are present in the raw data.
    """

    xrefs = [
        ('RefSeq', refseq_id),
        ('GtRNAdb', gtrnadb_id),
        ('Ensembl', ensembl_id),
    ]
    data = {}
    for name, func in xrefs:
        value = func(raw)
        if not value:
            continue
        data[name] = func(raw)
    return data


def as_partial_entry(raw):
    """
    Create a partial Entry. This will contain everything needed to build an
    Entry but a sequence.
    """
    return PartialEntry(
        primary_id=raw['symbol'],
        accession=raw['hgnc_id'],
        ncbi_tax_id=taxid,
        database='HGNC',
        exons=[],
        rna_type=rna_type(raw),
        url='',
        seq_version=1,
        note_data={},
        xref_data=xref_data(raw),
        species=species(raw),
        common_name=common_name(raw),
        lineage=lineage(raw),
        gene=raw['symbol'],
        organelle=organelle(raw),
        description=description(raw),
        gene_synonyms=gene_symbols(raw),
        references=references(raw),
    )
