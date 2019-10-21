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

import json
import pickle
import logging
import collections as coll

import attr

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import publications as pubs

from . import helpers

LOGGER = logging.getLogger(__name__)


def as_entry(row, reference_mapping):
    return data.Entry(
        primary_id=helpers.primary_id(row),
        accession=helpers.accession(row),
        ncbi_tax_id=helpers.taxid(row),
        database='PDBE',
        sequence=helpers.sequence(row),
        regions=[],
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


def as_mapping(report):
    mapping = coll.defaultdict(list)
    for row in report:
        ref = attr.asdict(helpers.as_reference(row))
        mapping[helpers.reference_mapping_id(row)].append(ref)
    return mapping


def as_reference_mapping(raw):
    mapping = coll.defaultdict(list)
    for pdb_id, refs in raw.items():
        for ref in refs:
            pub = attr.asdict(helpers.as_reference(ref))
            mapping[pdb_id.upper()].append(pub)
    return mapping


def as_descriptions(report):
    disqualified = {'mRNA': 0}

    descriptions = []
    for row in report:
        if helpers.is_mrna(row):
            disqualified['mRNA'] += 1

        if helpers.is_ncrna(row):
            descriptions.append(row)

    LOGGER.info('Disqualified %i mRNA chains', disqualified['mRNA'])
    return descriptions


def as_entries(raw_data, reference_mapping):
    for raw in raw_data:
        try:
            yield as_entry(raw, reference_mapping)
        except helpers.InvalidSequence as err:
            LOGGER.warn("Invalid sequence")
            LOGGER.exception(err)


def parse(handle, reference_handle):
    raw_data = json.load(handle)
    reference_mapping = pickle.load(reference_handle)
    return as_entries(raw_data, reference_mapping)
