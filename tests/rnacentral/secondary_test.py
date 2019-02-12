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

from rnacentral_pipeline.rnacentral import secondary as sec


@pytest.mark.parametrize('directory,count', [
    ('data/secondary', 1),
])
def test_can_process_a_directory(directory, count):
    assert len(list(sec.process_directory(directory))) == count


def test_can_produce_reasonable_data():
    val = list(sec.process_directory('data/secondary'))
    upi, model, dot, svg = val[0]
    assert upi == 'URS00000F9D45_9606'
    assert model == 'd.5.e.H.sapiens.2'
    assert len(dot) == 121
    assert svg.startswith('<svg')
    assert '\n' not in svg
