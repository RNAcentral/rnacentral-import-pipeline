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

import pytest

from rnacentral_pipeline.databases.psi_mi import data


@pytest.mark.parametrize('raw,ident', [
    ('psi-mi:"MI:0396"(predetermined participant)', data.Identifier('psi-mi', 'MI:0396', 'predetermined participant')),
    ('uniprotkb:O43426', data.Identifier('uniprotkb', 'O43426', None)),
    ('uniprotkb:P34708-1', data.Identifier('uniprotkb', 'P34708-1', None)),
    ('go:"GO:0005525"(GTP binding)', data.Identifier('go', "GO:0005525", 'GTP binding')),
    ('rcsb pdb:1DYN', data.Identifier('rcsb pdb', '1DYN', None)),
])
def test_can_bulid_ids(raw, ident):
    assert data.Identifier.build(raw) == ident

