# -*- coding: utf-8 -*-

"""
Copyright [2009-curren] EMBL-European Bioinformatics Institute
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

import csv
import operator as op

from rnacentral_pipeline.rnacentral import lookup

from .data import Context
from .data import KnownSequence

QUERY = """
SELECT
    pre.id as rna_id,
    pre.rna_type,
    COALESCE(rna.seq_short, rna.seq_long) as sequence,
    pre.description
from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
where
    pre.id in %s
"""


def parse_urs(handle, field: str):
    getter = op.itemgetter(field)
    reader = csv.DictReader(handle, delimiter='\t')
    return {getter(r) for r in reader}


def write(handle, db_url, field_name, output):
    data = parse_urs(handle, field_name)
    lookup.write_mapping(data, db_url, output, key='rna_id')


def load(handle):
    raw = lookup.load_mapping(handle)
    data = {}
    for key, value in raw.items():
        data[key] = KnownSequence(**value)
    return data
