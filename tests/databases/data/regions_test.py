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

from rnacentral_pipeline.databases.data.regions import *

@pytest.mark.parametrize('raw,expected', [
    (1, Strand.forward),
    (1.0, Strand.forward),
    ('1', Strand.forward),
    ('+', Strand.forward),
    (-1, Strand.reverse),
    ('-1', Strand.reverse),
    ('-', Strand.reverse),
    (0, Strand.unknown),
    (0.0, Strand.unknown),
    ('.', Strand.unknown),
])
def test_can_convert_a_strand(raw, expected):
    assert Strand.build(raw) is expected


@pytest.mark.parametrize('raw', [
    1.1,
    -0.2,
    'bob',
])
def test_fails_with_bad_strands(raw):
    with pytest.raises(UnknownStrand):
        Strand.build(raw)


@pytest.mark.parametrize('name,expected', [
    ('0-start, half-open', 
     CoordinateSystem(
         basis=CoordianteStart.zero, 
         close_status=CloseStatus.open)
     ),
    ('1-start, fully-closed', 
     CoordinateSystem(
         basis=CoordianteStart.one,
         close_status=CloseStatus.closed)
     ),
])
def test_can_build_coordinate_system_from_name(name, expected):
    assert CoordinateSystem.from_name(name) == expected
    assert CoordinateSystem.from_name(name).name() == name


@pytest.mark.parametrize('name,exon,expected', [
    ('0-start, half-open', Exon(start=10, stop=12), Exon(start=11, stop=11)),
    ('1-start, fully-closed', Exon(start=10, stop=12), Exon(start=10, stop=12)),
])
def test_can_correctly_normalize_an_exon(name, exon, expected):
    assert CoordinateSystem.from_name(name).normalize(exon) == expected


@pytest.mark.parametrize('strand,coordinate_system,expected', [
    (Strand.forward, CoordinateSystem.from_name('0-start, half-open'), '@4/11-21,31-39:+'),
    (Strand.reverse, CoordinateSystem.from_name('0-start, half-open'), '@4/11-21,31-39:-'),
    (Strand.forward, CoordinateSystem.from_name("1-start, fully-closed"), '@4/10-22,30-40:+'),
    (Strand.reverse, CoordinateSystem.from_name("1-start, fully-closed"), '@4/10-22,30-40:-'),
])
def test_can_generate_correct_region_names(strand, coordinate_system, expected):
    region = SequenceRegion(
        assembly_id='GRCh38',
        chromosome='4',
        strand=strand,
        exons=[
            Exon(start=10, stop=22),
            Exon(start=30, stop=40),
        ],
        coordinate_system=coordinate_system,
    )
    assert region.name() == expected


@pytest.mark.parametrize('strand,coordinate_system,expected', [
    (Strand.forward, CoordinateSystem.from_name('0-start, half-open'), [
        ['', '@4/31-39:+', '4', 1, 'GRCh38', 1, 31, 39]
    ]),
    (Strand.reverse, CoordinateSystem.from_name('0-start, half-open'), [
        ['', '@4/31-39:-', '4', -1, 'GRCh38', 1, 31, 39]
    ]),
    (Strand.forward, CoordinateSystem.from_name("1-start, fully-closed"), [
        ['', '@4/30-40:+', '4', 1, 'GRCh38', 1, 30, 40]
    ]),
    (Strand.reverse, CoordinateSystem.from_name("1-start, fully-closed"), [
        ['', '@4/30-40:-', '4', -1, 'GRCh38', 1, 30, 40]
    ]),
])
def test_can_generate_correct_writeable(strand, coordinate_system, expected):
    region = SequenceRegion(
        assembly_id='GRCh38',
        chromosome='4',
        strand=strand,
        exons=[Exon(start=30, stop=40)],
        coordinate_system=coordinate_system,
    )
    assert list(region.writeable('')) == expected
