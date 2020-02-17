# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.ncbi.gene import fetch
from rnacentral_pipeline.databases.ncbi.gene import helpers


@pytest.fixture(scope='module')
def simple_sequences():
    with open('data/ncbi_gene/simple.txt') as raw:
        ncrnas = helpers.ncrnas(raw)
        return fetch.sequences(ncrnas)


def test_can_get_expected_sequences(simple_sequences):
    assert len(simple_sequences) == 11


@pytest.mark.parametrize('gene_id,sequence_id', {
    ('108869573', 'NC_024802.1'),
    ('107303968', 'XR_001549991.1'),
})
def test_can_fetch_correct_sequence_data(simple_sequences, gene_id, sequence_id):
    assert simple_sequences[gene_id].id == sequence_id


@pytest.mark.skip(reason="Too long right now")
def test_can_fetch_ncrnas():
    # This was the count as of 2019-09-20, but I assume it will always increase
    # from here on out.
    total = sum(1 for ncrna in fetch.raw())
    assert total >= 2534341
