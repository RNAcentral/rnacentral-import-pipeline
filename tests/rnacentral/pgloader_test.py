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

from rnacentral_pipeline.rnacentral import pgloader


@pytest.mark.parametrize('output,expected', [
    ('data/pgloader/failed.txt', False),
    ('data/pgloader/success.txt', True),
])
def test_can_validate_output(output, expected):
    with open(output, 'r') as handle:
        assert pgloader.validate(handle) is expected
