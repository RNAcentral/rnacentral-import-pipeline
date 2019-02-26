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

from rnacentral_pipeline.databases.data import regions


@pytest.mark.parametrize('raw,expected', [
    (1, 1),
    (1.0, 1),
    ('1', 1),
    ('+', 1),
    (-1, -1),
    ('-1', -1),
    ('-', -1),
    (0, 0),
    (0.0, 0),
    ('.', 0),
])
def test_can_convert_a_strand(raw, expected):
    assert regions.as_strand(raw) == expected


@pytest.mark.parametrize('raw', [
    1.1,
    -0.2,
    'bob',
])
def test_fails_with_bad_strands(raw):
    with pytest.raises(regions.UnknownStrand):
        regions.as_strand(raw)
