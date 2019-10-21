# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from .. import helpers


@pytest.mark.parametrize('rna_id,expected', [
    ('URS000013F331_9606', {'Eukaryota'}),
    ('URS0000197DAA_77133', {'Bacteria'}),
    ('URS00001A89B6_137771', {'Eukaryota'}),
    ('URS00001FBB61_111789', set()),
    ('URS00002C6CD1_6239', {'Eukaryota'}),
    ('URS000040E1F7_562', {'Bacteria'}),
    ('URS00004761A1_660924', {'Eukaryota'}),
    ('URS0000759CF4_9606', {'Eukaryota'}),
    ('URS0000767631_155900', set()),
    ('URS00008089CA_224308', {'Bacteria'}),
    ('URS0000871A17_32630', set()),
    ('URS00008ABD77_155900', set()),
    ('URS0000D664B8_12908', set()),
])
def test_can_get_correct_domains(rna_id, expected):
    assert helpers.load_data(rna_id).domains() == expected


@pytest.mark.parametrize('rna_id,expected', [
    ('URS00002C6CD1_6239', False),
    ('URS000031188A_9606', True),
    ('URS000031415B_9606', True),
    ('URS00003EA5AC_9606', True),
    ('URS0000010837_7227', False),
    ('URS0000767631_155900', False),
])
def test_can_detect_if_mitochondrial(rna_id, expected):
    assert helpers.load_data(rna_id).is_mitochondrial() is expected


@pytest.mark.parametrize('rna_id,expected', [
    ('URS00002C6CD1_6239', False),
    ('URS00004DCBD3_3702', True),
    ('URS00006EF97C_3702', True),
    ('URS0000506D7B_100272', True),
    ('URS0000767631_155900', True),
])
def test_can_detect_if_chloroplast(rna_id, expected):
    assert helpers.load_data(rna_id).is_chloroplast() is expected


@pytest.mark.skip()
def test_can_correctly_load_hgnc_data():
    pass


@pytest.mark.skip()
def test_can_correctly_load_generic_data():
    pass
