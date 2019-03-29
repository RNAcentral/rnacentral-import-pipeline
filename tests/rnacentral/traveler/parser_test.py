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

from rnacentral_pipeline.rnacentral.traveler import parser as sec
from rnacentral_pipeline.databases.helpers.hashes import md5


@pytest.mark.parametrize('directory,count', [
    ('data/traveler/simple', 1),
])
def test_can_process_a_directory(directory, count):
    assert len(list(sec.models(directory))) == count


def test_can_produce_reasonable_data():
    val = list(sec.models('data/traveler/simple'))
    assert attr.asdict(val[0]) == attr.asdict(sec.TravelerResult(
        urs='URS00000F9D45_9606',
        model_id='d.5.e.H.sapiens.2',
        directory='data/traveler/simple',
        overlap_count=0,
        ribotyper=sec.ribotyper.Result(
            target='URS00000F9D45_9606',
            status='PASS',
            length=1588,
            fm=1,
            fam='SSU',
            domain='Bacteria',
            model='d.16.b.C.perfringens',
            strand=1,
            ht=1,
            tscore=1093.0,
            bscore=1093.0,
            bevalue=0.0,
            tcov=0.999,
            bcov=0.999,
            bfrom=3,
            bto=1588,
            mfrom=3,
            mto=1512,
        )))

def test_can_extract_expected_svg_data():
    val = list(sec.models('data/traveler/simple'))
    svg = val[0].svg()
    assert '\n' not in svg
    assert svg.startswith('<svg')
    assert md5(svg.encode()) == '2204b2f0ac616b8366a3b5f37aa123b8'


def test_can_extract_expected_dot_bracket_data():
    val = list(sec.models('data/traveler/simple'))
    assert len(val) == 1
    print(val[0])
    assert val[0].dot_bracket() == '(((((((((....((((((((.....((((((............))))..))....)))))).)).(((((......((.((.(((....))))).)).....))))).)))))))))...'
    assert val[0].basepair_count == 35
