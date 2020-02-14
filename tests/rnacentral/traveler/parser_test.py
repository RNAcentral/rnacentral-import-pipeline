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
    ('data/traveler/rfam', 2),
    ('data/traveler/ribovision', 2),
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

def test_produces_valid_data_for_rfam():
    val = list(sec.models('data/traveler/rfam'))
    assert attr.asdict(val[0]) == attr.asdict(sec.TravelerResult(
        urs='URS0000A7635A',
        model_id='RF00162',
        directory='data/traveler/rfam',
        overlap_count=0,
        ribotyper=None,
        is_rfam=True,
    ))

    writeable = val[0].writeable() 

    # This is too long to include in the file so I just test the hash and delete
    # it
    assert md5(writeable[3].encode()) == '9504c4b9a1cea77fa2c4ef8082d7b996'
    del writeable[3]
    
    assert writeable == [
       'URS0000A7635A' ,
        'RF00162',
        '((((((((......(((...(((.....)))......))).(((.(((......))))))........((((......))))...)))))))).',
        0,
        24,
        None,
        None,
        None,
        None,
        None,
    ]


def test_parses_ribovision_results():
    vals = list(sec.models('data/traveler/ribovision'))
    val = next(v for v in vals if v.urs == 'URS0000C5FF65')
    val = attr.asdict(val)
    from pprint import pprint
    pprint(val)
    assert val == attr.asdict(sec.TravelerResult(
        urs='URS0000C5FF65',
        model_id='EC_LSU_3D',
        directory='data/traveler/ribovision',
        overlap_count=3,
        ribotyper=sec.ribotyper.Result(
            target='URS0000C5FF65',
            status='PASS',
            length=2892,
            fm=1,
            fam='LSU',
            domain='Bacteria',
            model='EC_LSU_3D',
            strand=1,
            ht=1,
            tscore=2197.9,
            bscore=2197.9,
            bevalue=0.0,
            tcov=0.999,
            bcov=0.999,
            bfrom=2,
            bto=2890,
            mfrom=2,
            mto=2902,
        )))


@pytest.mark.parametrize('directory,index,md5_hash', [
    ('data/traveler/simple', 0, '2204b2f0ac616b8366a3b5f37aa123b8'),
    ('data/traveler/rfam', 0, '9504c4b9a1cea77fa2c4ef8082d7b996'),
])
def test_can_extract_expected_svg_data(directory, index, md5_hash):
    val = list(sec.models(directory))
    svg = val[index].svg()
    assert '\n' not in svg
    assert svg.startswith('<svg')
    assert md5(svg.encode()) == md5_hash


@pytest.mark.parametrize('directory,index,secondary,bp_count', [
    ('data/traveler/simple', 0, '(((((((((....((((((((.....((((((............))))..))....)))))).)).(((((......((.((.(((....))))).)).....))))).)))))))))...', 35),
    ('data/traveler/rfam', 1, '((((((((......(((...(((.....)))......))).(((.(((......))))))........((((......))))...)))))))).', 24),
])
def test_can_extract_expected_dot_bracket_data(directory, index, secondary,
                                               bp_count):
    val = list(sec.models(directory))
    assert val[index].dot_bracket() == secondary
    assert val[index].basepair_count == bp_count
