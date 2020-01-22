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
from pathlib import Path

from rnacentral_pipeline.rnacentral.traveler import parser
from rnacentral_pipeline.rnacentral.traveler import data

from rnacentral_pipeline.databases.helpers.hashes import md5


@pytest.mark.parametrize('directory,source,count', [
    ('data/traveler/crw', data.Source.crw, 1),
    ('data/traveler/rfam', data.Source.rfam, 2),
    ('data/traveler/ribovision', data.Source.ribovision, 2),
])
def test_can_process_a_directory(directory, source, count):
    path = Path(directory)
    assert len(list(parser.parse(source, path))) == count


def test_can_produce_reasonable_data():
    val = list(parser.parse(data.Source.crw, Path('data/traveler/crw')))
    assert attr.asdict(val[0]) == attr.asdict(data.TravelerResult(
        urs='URS00000F9D45_9606',
        model_id='d.5.e.H.sapiens.2',
        basepath=Path('data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2'),
        source=data.Source.crw,
        ribovore=data.RibovoreResult(
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


@pytest.mark.parametrize('directory,source,index,attr,expected', [
    ('data/traveler/crw', data.Source.crw, 0, 'dot_bracket_path', 'data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.fasta'),
    ('data/traveler/crw', data.Source.crw, 0, 'svg_path', 'data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg'),
])
def test_can_produce_correct_paths(directory, source, index, attr, expected):
    val = list(parser.parse(source, Path(directory)))
    assert getattr(val[index], attr) == Path(expected)


def test_gets_correct_count():
    val = list(parser.parse(data.Source.rfam, Path('data/traveler/rfam')))
    v = next(v for v in val if v.urs == 'URS0000A7635A')
    assert v.overlap_count() == 0
    assert v.basepair_count() == 24


def test_produces_valid_data_for_rfam():
    val = parser.parse(data.Source.rfam, Path('data/traveler/rfam'))
    v = next(v for v in val if v.urs == 'URS0000A7635A')
    assert attr.asdict(v) == attr.asdict(data.TravelerResult(
        urs='URS0000A7635A',
        model_id='RF00162',
        basepath=Path('data/traveler/rfam/RF00162/URS0000A7635A'),
        source=data.Source.rfam,
        ribovore=None,
    ))

    assert v.svg_path == Path('data/traveler/rfam/RF00162/URS0000A7635A.colored.svg')
    writeable = v.writeable()

    # This is too long to include in the file so I just compare to the file it
    # is meant to read.
    with Path('data/traveler/rfam/RF00162/URS0000A7635A.colored.svg').open('r') as raw:
        d = raw.read().replace('\n', '')
        assert d == writeable[3]
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
        '',
    ]


def test_parses_ribovision_results():
    vals = parser.parse(data.Source.ribovision, Path('data/traveler/ribovision'))
    val = next(v for v in vals if v.urs == 'URS0000C5FF65')
    assert attr.asdict(val) == attr.asdict(data.TravelerResult(
        urs='URS0000C5FF65',
        model_id='EC_LSU_3D',
        basepath=Path('data/traveler/ribovision/URS0000C5FF65-EC_LSU_3D'),
        source=data.Source.ribovision,
        ribovore=data.RibovoreResult(
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


@pytest.mark.parametrize('directory,source,urs,md5_hash', [
    ('data/traveler/crw', data.Source.crw, 'URS00000F9D45_9606', '2204b2f0ac616b8366a3b5f37aa123b8'),
    ('data/traveler/rfam', data.Source.rfam, 'URS0000A7635A', '9504c4b9a1cea77fa2c4ef8082d7b996'),
])
def test_can_extract_expected_svg_data(directory, source, urs, md5_hash):
    val = list(parser.parse(source, Path(directory)))
    svg = next(v for v in val if v.urs == urs).svg()
    assert '\n' not in svg
    assert svg.startswith('<svg')
    assert md5(svg.encode()) == md5_hash


@pytest.mark.parametrize('directory,source,urs,secondary,bp_count', [
    ('data/traveler/crw', data.Source.crw, 'URS00000F9D45_9606', '(((((((((....((((((((.....((((((............))))..))....)))))).)).(((((......((.((.(((....))))).)).....))))).)))))))))...', 35),
    ('data/traveler/rfam', data.Source.rfam, 'URS0000A7635B', '((((((((......(((...(((.....)))......))).(((.(((......))))))........((((......))))...)))))))).', 24),
])
def test_can_extract_expected_dot_bracket_data(directory, source, urs, secondary,
                                               bp_count):
    val = parser.parse(source, Path(directory))
    v = next(v for v in val if v.urs == urs)
    assert v.dot_bracket() == secondary
    assert v.basepair_count() == bp_count
