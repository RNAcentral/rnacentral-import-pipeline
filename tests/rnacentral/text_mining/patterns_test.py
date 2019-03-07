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

from rnacentral_pipeline.rnacentral.text_mining import patterns


@pytest.mark.parametrize('pats,filename,found', [
    ([r'HO[TX]', r'rs\d+'], 'data/text-mining/PMID:29141248.txt', 
     {'HOX', 'rs920778', 'rs4759314', 'rs1899663'}),
    (['XIST'], 'data/text-mining/PMID:30541551.txt', set()),
])
def test_can_find_expected_matches_in_text(pats, filename, found):
    with open(filename, 'r') as raw:
        assert set(patterns.matches(pats, raw)) == found
