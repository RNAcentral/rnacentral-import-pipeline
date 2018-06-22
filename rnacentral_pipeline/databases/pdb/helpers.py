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
import csv
import logging
import itertools as it
import collections as coll

import requests

from databases.helpers import phylogeny as phy
from databases.data import Reference

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

ALLOWED = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


class InvalidSequence(Exception):
    pass


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return it.izip_longest(*args, fillvalue=fillvalue)


def rna_containing_pdb_ids():
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """
    query = """
    <orgPdbQuery>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <containsProtein>I</containsProtein>
    <containsDna>I</containsDna>
    <containsRna>Y</containsRna>
    <containsHybrid>I</containsHybrid>
    </orgPdbQuery>
    """
    url = 'http://www.rcsb.org/pdb/rest/search'

    # without this header the request is redirected incorrectly
    request = requests.post(
        url,
        data=query,
        headers={'content-type': 'application/x-www-form-urlencoded'}
    )

    if request.status_code == 200:
        pdb_ids = request.text.rstrip().split('\n')
    else:
        pdb_ids = None

    return pdb_ids


def custom_report(pdb_ids, fields):
    """
    Get custom report about PDB files in a tabular format.
    """

    url = 'http://www.rcsb.org/pdb/rest/customReport.csv'
    data = {
        'pdbids': ','.join(pdb_ids),
        'customReportColumns': ','.join(fields),
        'format': 'csv',
        'service': 'wsfile',  # Use actual CSV files
    }

    response = requests.post(url, data=data)
    response.raise_for_status()
    lines = response.text.split('\n')
    return csv.DictReader(lines, delimiter=',', quotechar='"')


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
    parts = row['citationAuthor'].split(',')
    grouped = grouper(parts, 2, '')
    authors = ','.join(''.join(p) for p in grouped)

    pmid = row['pubmedId']
    if pmid:
        pmid = int(pmid)
    else:
        pmid = None

    return Reference(
        authors=authors,
        location=row['journalName'],
        title=row['title'],
        pmid=pmid,
        doi=row['doi'],
    )


def reference_mapping(pdb_ids):
    """
    Get literature citations for each PDB file.
    Return a dictionary that looks like this:
    {
        '1S72': {'structureId': '1S72', etc}
    }
    """

    report = custom_report(pdb_ids, [
        'structureId',
        'citationAuthor',
        'firstPage',
        'lastPage',
        'journalName',
        'title',
        'volumeId',
        'publicationYear',
        'pubmedId',
        'pmc',
        'doi',
    ])

    mapping = coll.defaultdict(list)
    for row in report:
        mapping[reference_mapping_id(row)].append(as_reference(row))
    return mapping


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
