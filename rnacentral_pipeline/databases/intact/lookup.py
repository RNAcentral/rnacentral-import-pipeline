# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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

import pickle
import operator as op

from rnacentral_pipeline import utils
from rnacentral_pipeline.rnacentral import lookup


QUERY = """
select
    pre.id as id,
    pre.rna_type,
    COALESCE(rna.seq_short, rna.seq_long) as sequence,
    pre.description
from rnc_rna_precomputed pre
join rna on rna.upi = pre.upi
where
    pre.id in %s
"""


def ids(interactions):
    getter = op.attrgetter("urs_taxid")
    return {getter(r) for r in interactions}


def mapping(db_url, data):
    mapping = lookup.as_mapping(db_url, ids(data), QUERY)
    for value in mapping.values():
        value["sequence"] = value["sequence"].replace("U", "T")
    return mapping
