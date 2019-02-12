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

import attr
import pytest

from rnacentral_pipeline.databases import pfam

@pytest.fixture
def simple():
    with open('data/pfam/pfam.tsv') as raw:
        yield list(pfam.parse(raw))


def test_can_parse_very_simple_data(simple):
    assert len(simple) == 1
    assert attr.asdict(simple[0]) == attr.asdict(pfam.Hit(
        seq=pfam.SequenceComponent(
            urs='ENA|Z11840.1:2823..6850:misc_RNA|Z11840.1:2823..6850:misc_RNA.1:0',
            start=173,
            stop=274,
            strand=1,
        ),
        model=pfam.ModelComponent(
            model_id='PF01085.18',
            name='HH_signal',
            clan_id='CL0170',
            start=1,
            stop=95,
            length=161,
        ),
        bit_score=95.8,
        e_value=2e-27,
    ))
