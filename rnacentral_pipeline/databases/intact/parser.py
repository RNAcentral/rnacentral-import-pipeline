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

import operator as op
import itertools as it

from rnacentral_pipeline.databases.psi_mi import tab

from . import lookup
from . import helpers

IGNORE = {
    'URS00001346E7_559292',
    'URS00002C74E0_559292',
    'URS00002C74E0 _559292',
    'URS0000459091_559191',
    'URS000057FD8B-559292',  # Is known, but formatted wrong
    'URS0000643233_559292',
    'URS00008FED3E_559292',
}


def parse_interactions(handle):
    data = tab.parse(handle)
    data = filter(op.methodcaller('involves_rnacentral'), data)
    data = filter(lambda i: i.urs_taxid.startswith('URS'), data)
    data = filter(lambda i: i.urs_taxid not in IGNORE, data)
    data = filter(lambda i: "_" in i.urs_taxid, data)
    return data


def parse(handle, db_url):
    key = op.attrgetter('urs_taxid')
    interactions = sorted(parse_interactions(handle), key=key)
    mapping = lookup.mapping(db_url, interactions)
    interactions = sorted(interactions, key=key)
    grouped = it.groupby(interactions, key)
    for urs_taxid, interactions in grouped:
        if urs_taxid not in mapping:
            raise ValueError("Found no sequence info for %s" % urs_taxid)

        info = mapping[urs_taxid]
        entry = helpers.as_entry(urs_taxid, list(interactions), info)
        yield entry
