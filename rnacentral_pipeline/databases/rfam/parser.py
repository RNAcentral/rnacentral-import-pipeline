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

import itertools as it
import json

from ..data import Entry
from . import helpers, utils


def as_entry(data):
    """
    Turn an entry from the JSON file into a Entry object for writing.
    """

    location_range = helpers.feature_location_endpoints(data)
    return Entry(
        primary_id=helpers.primary_id(data),
        accession=helpers.accession(data),
        ncbi_tax_id=helpers.taxid(data),
        database="RFAM",
        sequence=helpers.sequence(data),
        regions=[],
        rna_type=helpers.rna_type(data),
        url=helpers.url(data),
        note_data=helpers.note(data),
        species=helpers.species(data),
        common_name=helpers.common_name(data),
        lineage=helpers.lineage(data),
        optional_id=helpers.optional_id(data),
        product=helpers.product(data),
        parent_accession=helpers.parent_accession(data),
        project="RFAM",
        experiment=helpers.experiment(data),
        description=helpers.description(data),
        mol_type=helpers.mol_type(data),
        seq_version=helpers.seq_version(data),
        is_composite="N",
        location_start=location_range[0],
        location_end=location_range[1],
        references=helpers.references(data),
    )


def parse(handle, mapping_file):
    """
    Parse the JSON file of Rfam data and produce a generator of all Entry
    objects in the file.
    """

    data = json.load(handle)
    return map(lambda e: as_entry(e), data)
