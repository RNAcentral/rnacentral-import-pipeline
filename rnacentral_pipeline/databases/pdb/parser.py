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
import typing as ty
import collections as coll

import attr

from rnacentral_pipeline.databases import data

from rnacentral_pipeline.databases.pdb.data import ChainInfo
from rnacentral_pipeline.databases.pdb import helpers

LOGGER = logging.getLogger(__name__)


def as_entry(info: ChainInfo, reference_mapping):
    return data.Entry(
        primary_id=info.pdb_id,
        accession=info.accession,
        ncbi_tax_id=helpers.taxid(info),
        database='PDBE',
        sequence=helpers.sequence(info),
        regions=[],
        rna_type=helpers.rna_type(row),
        url=helpers.url(info),
        seq_version='1',
        note_data=helpers.note_data(row),
        xref_data=helpers.xref_data(row),
        product=helpers.product(row),
        optional_id=info.chain_id,
        description=helpers.description(info),
        species=helpers.species(info),
        lineage=helpers.lineage(info),
        parent_accession=info.pdb_id,
        references=helpers.references_for(info, reference_mapping),
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


def as_entries(raw_data: ty.List[ChainInfo], reference_mapping) -> ty.Iterator[data.Entry]:
    disqualified = {'mRNA': 0, "other": 0}
    for raw in raw_data:
        if helpers.is_mrna(raw):
            LOGGER.debug("Disqualifing %s", raw)
            disqualified['mRNA'] += 1
            continue

        if not helpers.is_ncrna(raw):
            LOGGER.debug("Skipping %s", raw)
            disqualified["other"] += 1
            continue

        try:
            yield as_entry(raw, reference_mapping)
        except helpers.InvalidSequence as err:
            LOGGER.warn("Invalid sequence")
            LOGGER.exception(err)
    LOGGER.info('Disqualified %i mRNA chains', disqualified['mRNA'])


def parse(handle, reference_handle):
    raw_data = json.load(handle)
    reference_mapping = pickle.load(reference_handle)
    return as_entries(raw_data, reference_mapping)
