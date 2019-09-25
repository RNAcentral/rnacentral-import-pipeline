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


def test_can_get_expected_sequences():
    with open('data/ncbi_gene/simple.txt') as raw:
        ncrnas = list(helpers.ncrnas(raw))
        seqs = fetch.sequences(ncrnas)
        assert len(seqs) == len(ncrnas)
        assert set(seqs.keys()) == set(n['GeneID'] for n in ncrnas)


def test_can_fetch_data():
    assert fetch.raw().read()


def test_can_fetch_ncrnas():
    data = fetch.ncrnas()
    # This was the count as of 2019-09-20, but I assume it will always increase
    # from here on out.
    assert sum(1 for n in fetch.ncrnas()) >= 2534341
