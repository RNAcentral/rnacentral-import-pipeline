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

from . import utils
from databases import data

LOGGER = logging.getLogger(__name__)


def chain_descriptions(pdb_ids):
    """
    Get per-chain information about each RNA sequence.
    Return a dictionary that looks like this:
    {
        '1S72_A': {'structureId': '1S72', 'chainId': '0' etc}
    }
    """

    report = utils.custom_report(pdb_ids, [
        'structureId',
        'chainId',
        'structureTitle',
        'experimentalTechnique',
        'releaseDate',
        'ndbId',
        'emdbId',
        'classification',
        'entityId',
        'sequence',
        'chainLength',
        'db_id',
        'db_name',
        'entityMacromoleculeType',
        'source',
        'taxonomyId',
        'compound',
        'resolution',
    ])

    disqualified = {'mRNA': 0}

    for row in report:
        # skip proteins
        if 'RNA' not in row['entityMacromoleculeType']:
            continue

        # skip mRNAs
        if re.search('mRNA', row['compound'], re.IGNORECASE) and \
           not re.search('tmRNA', row['compound'], re.IGNORECASE):
            disqualified['mRNA'] += 1
            continue

        yield row

    LOGGER.info('Disqualified %i mRNA chains', disqualified['mRNA'])


def as_entry(row, reference_mapping):
    return data.Entry(
        primary_id=utils.primary_id(row),
        accession=utils.accession(row),
        ncbi_tax_id=utils.taxid(row),
        database='PDBE',
        sequence=utils.sequence(row),
        exons=[],
        rna_type=utils.rna_type(row),
        url=utils.url(row),
        seq_version='1',
        note_data=utils.note_data(row),
        xref_data=utils.xref_data(row),
        product=utils.product(row),
        optional_id=utils.optional_id(row),
        description=utils.description(row),
        parent_accession=utils.parent_accession(row),
        references=utils.references_for(row, reference_mapping),
        location_start=utils.location_start(row),
        location_end=utils.location_end(row),
    )


def as_entries(pdb_ids):
    reference_mapping = utils.reference_mapping(pdb_ids)
    for result in chain_descriptions(pdb_ids):
        yield as_entry(result, reference_mapping)
