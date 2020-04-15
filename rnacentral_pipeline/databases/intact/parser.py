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

import csv
import operator as op
import itertools as it

from rnacentral_pipeline.writers import write_entries
from rnacentral_pipeline.databases.psi_mi import tab

from . import lookup


def parse(handle, db_url):
    key = op.attrgetter('rnacentral_iteractor')
    interactions = tab.parse(handle)
    interactions = filter(op.methodcaller('involves_rnacentral'), rows)
    interactions = sorted(interactions, key=key)
    mapping = lookup.mapping(db_url, rows)
    grouped = it.groupby(interactions, key)
    for urs_taxid, interactions in grouped:
        interactions = list(interactions)
        entry = helpers.as_entry(urs_taxid, interactions, mapping)
        yield entry
