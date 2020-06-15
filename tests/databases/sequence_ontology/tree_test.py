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

from rnacentral_pipeline.databases.sequence_ontology import tree as so

import pytest


@pytest.mark.parametrize("so_term_id,expected", [
    ('rRNA', [('SO:0000655', 'ncRNA'), ('SO:0000252', 'rRNA')]),
    ('rRNA_18S', [
        ('SO:0000655', 'ncRNA'), 
        ('SO:0000252', 'rRNA'),
        ('SO:0000650', 'small_subunit_rRNA'),
        ('SO:0000407', 'rRNA_18S'),
    ]),
    ('lincRNA', [
        ('SO:0000655', 'ncRNA'), 
        ('SO:0001877', 'lnc_RNA'),
        ('SO:0001463', 'lincRNA'),
    ]),
])
def test_can_compute_some_simple_paths(so_term_id, expected):
    data = so.load_ontology(so.REMOTE_ONTOLOGY)
    assert so.rna_type_tree(data, so_term_id) == expected
