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

from rnacentral_pipeline.databases.ncbi.gene import helpers

@pytest.mark.paramterize('filename,count', [
    ('data/ncbi_gene/simple.txt', 3),
    ('data/ncbi_gene/ncrna.txt', 5),
])
def test_can_find_expected_ncrnas(filename, count):
    with open(filename) as raw:
        with len(list(helpers.ncrnas(raw))) == count


@pytest.mark.paramterize('filename', [
    ('data/ncbi_gene/proteins.txt'),
    ('data/ncbi_gene/empty.txt'),
])
def test_complains_if_no_ncrnas(filename):
    with pytest.raises(Exception):
        with open(filename) as raw:
            with list(helpers.ncrnas(raw))
