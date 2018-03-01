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

from . import helpers
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

    report = helpers.custom_report(pdb_ids, [
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
        primary_id=helpers.primary_id(row),
        accession=helpers.accession(row),
        ncbi_tax_id=helpers.taxid(row),
        database='PDBE',
        sequence=helpers.sequence(row),
        exons=[],
        rna_type=helpers.rna_type(row),
        url=helpers.url(row),
        seq_version='1',
        note_data=helpers.note_data(row),
        xref_data=helpers.xref_data(row),
        product=helpers.product(row),
        optional_id=helpers.optional_id(row),
        description=helpers.description(row),
        species=helpers.species(row),
        lineage=helpers.lineage(row),
        parent_accession=helpers.parent_accession(row),
        references=helpers.references_for(row, reference_mapping),
        location_start=helpers.location_start(row),
        location_end=helpers.location_end(row),
    )


def as_entries(pdb_ids):
    reference_mapping = helpers.reference_mapping(pdb_ids)
    for result in chain_descriptions(pdb_ids):
        yield as_entry(result, reference_mapping)
