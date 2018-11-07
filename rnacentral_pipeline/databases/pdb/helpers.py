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

import re
import logging
import itertools as it
import collections as coll

from rnacentral_pipeline.databases.helpers import phylogeny as phy
from rnacentral_pipeline.databases.helpers import publications as pubs
from rnacentral_pipeline.databases.data import Reference

RIBOSOMES = set([
    '5S',
    '5.8S',
    '16S',
    '18S',
    '23S',
    '28S',
    '30S',
    '40S',
    '60S',
    '80S',
])

URL = 'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}'

LOGGER = logging.getLogger(__name__)

ALLOWED = re.compile('^[ABCDFGHIKMNRSTVWXYU]+$', re.IGNORECASE)


class InvalidSequence(Exception):
    pass


def is_mrna(row):
    if re.search('mRNA', row['compound'], re.IGNORECASE) and \
       not re.search('tmRNA', row['compound'], re.IGNORECASE):
        return True
    return False


def is_ncrna(row):
    return 'RNA' in row['entityMacromoleculeType'] and not is_mrna(row)


def accession(row):
    """
    Generates and accession for the given entry. The accession is built from
    the structureId, chainId, and entityId.
    """

    # use entityId to ensure that the id is unique when chainIds
    # are only different in case ('A' and 'a')
    return '{structureId}_{chainId}_{entityId}'.format(
        structureId=row['structureId'],
        chainId=row['chainId'],
        entityId=row['entityId'],
    )


def sequence(row):
    """
    Fetches the sequence of the row as DNA.
    """
    sequence = row['sequence'].replace('U', 'T')
    # In many tRNA's there is a single last amino acid, so we ignore that and
    # check before excluding sequences.
    if not re.match(ALLOWED, sequence):
        raise InvalidSequence("%s appears to be mislabelled protein" % row)
    return sequence


def taxid(row):
    """
    Fetch the taxid from the row. This will deal with , or empty taxids.
    """

    # if no taxonomy id, use that of the synthetic construct
    if row['taxonomyId'] == '':
        return 32630  # synthetic construct

    # If there is a ',' in then that is synthetic
    if ',' in row['taxonomyId']:
        return 32630

    return int(row['taxonomyId'])


def as_reference(row):
    # This is pretty dirty but it should work assuming that each name
    # always has a ',' after both the first and last name.

    pmid = None
    if row['pubmed_id']:
        return pubs.reference(int(row['pubmed_id']))

    if row['doi']:
        doi = row['doi']
        if not doi.lower().startswith('doi:'):
            doi = 'doi:' + doi
        return pubs.reference(doi)

    authors = []
    for author in row['author_list']:
        if author['full_name']:
            authors.append(author['full_name'].replace(',', ''))
            continue

        current = []
        if author['last_name']:
            current.append(author['last_name'])
        if author['first_name']:
            current.append(author['first_name'])
        authors.append(' '.join(current))

    journal = row['journal_info']['pdb_abbreviation']

    return Reference(
        authors=', '.join(authors),
        location=journal,
        title=row['title'],
        pmid=pmid,
        doi=row['doi'],
    )


def references_for(row, mapping):
    return mapping[reference_mapping_id(row)]


def primary_id(row):
    return row['structureId']


def rna_type(row):
    compound = row['compound'].upper()
    for simple_type in ['tRNA', 'tmRNA', 'snRNA']:
        if simple_type.upper() in compound:
            return simple_type

    # rRNA
    for ribo_name in RIBOSOMES:
        if ribo_name in compound:
            return 'rRNA'

    # SRP
    if 'SRP' in compound:
        return 'SRP_RNA'

    # Ribozyme
    if 'RIBOZYME' in compound and 'HAMMERHEAD' not in compound:
        return 'ribozyme'

    # Hammerhead ribozyme
    if 'RIBOZYME' in compound and 'HAMMERHEAD' in compound:
        return 'hammerhead_ribozyme'

    # snoRNA
    if 'SNORNA' in compound:
        return 'ncRNA'

    return 'misc_RNA'


def url(row):
    """
    Generate a URL for a given result. It will point to the page for the whole
    structure.
    """
    return URL.format(pdb_id=row['structureId'].lower())


def xref_data(row):
    """
    Put NDB and EMDB xrefs in the db_xref field.
    """

    xref = coll.defaultdict(list)
    if row.get('ndbId', None):
        xref['NDB'].append(row['ndbId'])
    if row.get('db_name', None):
        db_name = row['db_name']
        if db_name != 'PDB':
            xref[db_name].append(row['db_id'])
    return dict(xref)


def note_data(row):
    fields = [
        'structureTitle',
        'experimentalTechnique',
        'resolution',
        'releaseDate',
    ]
    notes = {}
    for field in fields:
        if field in row and row[field]:
            notes[field] = row[field]
    return notes


def description(row, max_length=80):
    compound = row['compound'][:max_length] + \
                (row['compound'][max_length:] and '...')
    return '{compound} from {source} (PDB {pdb}, chain {chain})'.format(
        compound=compound,
        source=row['source'],
        pdb=row['structureId'],
        chain=row['chainId'],
    )


def product(row):
    return row['compound']


def optional_id(row):
    return row['chainId']


def reference_mapping_id(row):
    return row['structureId']


def parent_accession(row):
    return row['structureId']


def location_start(_):
    return 1


def location_end(row):
    return int(row['chainLength'])


def lineage(row):
    return phy.lineage(taxid(row))


def species(row):
    return phy.species(taxid(row))
