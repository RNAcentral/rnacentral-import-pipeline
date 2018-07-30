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
import logging

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases import mapping

LOGGER = logging.getLogger(__name__)

URL = 'https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id={id}'


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

REFSEQ_QUERY = """
select distinct
    acc.external_id as id,
    COALESCE(rna.seq_short, rna.seq_long) as sequence
from xref
join rnc_accessions acc on acc.accession = xref.ac
join rna on rna.upi = xref.upi
where
    xref.deleted = 'N'
    and xref.dbid = 9
    and xref.taxid = 9606
"""

ENSEMBL_QUERY = '''
select distinct
    split_part(acc.optional_id, '.', 1) as id,
    COALESCE(rna.seq_short, rna.seq_long) as sequence
from xref
join rnc_accessions acc on acc.accession = xref.ac
join rna on rna.upi = xref.upi
where
    xref.deleted = 'N'
    and xref.dbid = 25
    and xref.taxid = 9606
'''

GTRNADB_QUERY = """
select distinct
    acc.optional_id as id,
    COALESCE(rna.seq_short, rna.seq_long) as sequence
from xref
join rnc_accessions acc on acc.accession = xref.ac
join rna on rna.upi = xref.upi
where
    xref.deleted = 'N'
    and xref.dbid = 8
    and xref.taxid = 9606
"""


def known_refseq(psql_wrapper):
    """
    Get a Matcher that will matcher against all RefSeq human data.
    """
    return mapping.Matcher.from_iterable(
        'RefSeq',
        psql_wrapper.copy_to_iterable(REFSEQ_QUERY),
    )


def known_gtrnadb(psql_wrapper):
    """
    Get a Matcher that will match against all GtRNAdb human sequences.
    """
    return mapping.Matcher.from_iterable(
        'GtRNAdb',
        psql_wrapper.copy_to_iterable(GTRNADB_QUERY),
    )


def known_ensembl(psql_wrapper):
    """
    Get a Matcher that will match against known Ensembl data for humans.
    """
    return mapping.Matcher.from_iterable(
        'Ensembl',
        psql_wrapper.copy_to_iterable(ENSEMBL_QUERY),
    )


def known(dbconf):
    """
    Get a multi matcher for all databases we can match RefSeq data against.
    """

    psql_wrapper = psql.PsqlWrapper(dbconf)
    return mapping.MultiMatcher.build(
        known_refseq(psql_wrapper),
        known_gtrnadb(psql_wrapper),
        known_ensembl(psql_wrapper),
    )


def gene_symbols(raw):
    """
    Get a list of previous symbols this gene has been known as.
    """
    return raw.get('prev_symbol', [])


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
    return 'Homo sapiens (human) ' + raw['name'] + gene_description


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

    refs = []
    for pmid in raw.get('pubmed_id', []):
        try:
            refs.append(pubs.reference(pmid))
        except pubs.FailedPublicationId:
            LOGGER.warn("Could not fetch publications for %i", pmid)
        except pubs.UnknownPmid:
            LOGGER.error("Publication %i is not valid", pmid)
    return refs


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
    """
    This will generate an GtRNAdb ID given raw data from the
    """

    prefix = None
    patterns = {
        'tRNA': r'^TR(\S+)-(\S{3})(\d+-\d+)',
        'nmt-tRNA': r'^NMTR(\S+)-(\S{3})(\d+-\d+)',
    }
    for prefix, pattern in patterns.items():
        match = re.match(pattern, raw['symbol'])
        if match:
            break
    else:
        return None

    return '{prefix}-{aa}-{anticodon}-{range}'.format(
        prefix=prefix,
        aa=ONE_TO_THREE[match.group(1)],
        anticodon=match.group(2),
        range=match.group(3),
    )


def xref_data(raw):
    """
    Compute xrefs from HGNC. This will only contain the xrefs that we can map
    and are present in the raw data.
    """

    xrefs = [
        ('RefSeq', 'refseq_accession'),
        ('UCSC', 'ucsc_id'),
        ('Ensembl', 'ensembl_gene_id'),
        ('ENA', 'ena'),
        ('LNCipedia', 'lncipedia'),
    ]
    data = {}
    for name, item in xrefs:
        if item not in raw:
            continue
        value = raw[item]
        if not isinstance(value, list):
            value = [value]
        data[name] = value

    trna_id = gtrnadb_id(raw)
    if trna_id:
        data['GtRNAdb'] = [trna_id]

    return data


def as_partial_entry(raw):
    """
    Create a partial Entry. This will contain everything needed to build an
    Entry but a sequence.
    """
    return mapping.SequencelessEntry(
        primary_id=raw['symbol'],
        accession=raw['hgnc_id'],
        ncbi_tax_id=taxid(raw),
        database='HGNC',
        exons=[],
        rna_type=rna_type(raw),
        url=URL.format(id=raw['hgnc_id']),
        seq_version='1',
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
