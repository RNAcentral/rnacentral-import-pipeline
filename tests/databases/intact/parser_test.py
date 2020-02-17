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

from rnacentral_pipeline.databases.intact import parser


@pytest.mark.parametrize('filename,count', [
    ('data/intact/intact.txt', 12),
])
def test_can_parse_all_data(filename, count):
    with open(filename, 'r') as raw:
        data = list(parser.parse(raw))
        assert len(data) == count


@pytest.mark.skip()
def test_produces_correct_data():
    assert False, "NYI"
