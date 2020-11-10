# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from functools import lru_cache

import pytest

from rnacentral_pipeline.rnacentral.genes import build, data


@lru_cache(1000)
def load_data(assembly):
    return []


@lru_cache(1000)
def genes(assembly):
    locations = load_data(assembly)
    return genes.build(locations)


def genes_for(assembly, urs_taxid):
    found = genes(assembly)
    return [g for g in found if g.urs_taxid == urs_taxid]


@pytest.mark.parameterize('location_id,status', [
    (33989329, data.DataType.rejected),
    (33919630, data.DataType.rejected),
    (34007230, data.DataType.rejected),
    (34007229, data.DataType.rejected),
])
def test_selects_correct_representative(assembly, location_id, status):
    assert genes_for(assembly, location_id)
    # genes = genes_for(assembly, urs_taxid)
    # assert any(g.representative for g in genes) == representative
