# -*- coding: utf-8 -*-

# pylint: disable=no-member,missing-docstring,invalid-name

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

from rnacentral_pipeline.databases.data import regions
from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data

from .helpers import fetch_coord


def fetch_one(rna_id, assembly):
    coords = list(fetch_coord(rna_id, assembly))
    assert len(coords) == 1
    return coords[0]


@pytest.mark.parametrize('rna_id,assembly,count', [
    ('URS00008B37EC_9606', 'GRCh38', 1),
    ('URS00008C1914_9606', 'GRCh38', 1),
    ('URS00006683B1_281687', 'GRCh38', 0),
])
def test_can_fetch_all_coordinates_for_upi_taxid(rna_id, assembly, count):
    assert len(list(fetch_coord(rna_id, assembly))) == count


def test_can_find_correct_for_something_with_several_exons():
    found = fetch_one('URS00009BF201_9606', 'GRCh38')
    assert attr.asdict(found) == attr.asdict(data.Region(
        region_id='URS00009BF201_9606.0',
        rna_id='URS00009BF201_9606',
        region=regions.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='16',
            strand=-1,
            exons=(
                regions.Exon(start=14085, stop=14511),
                regions.Exon(start=14652, stop=14720),
                regions.Exon(start=15481, stop=15633),
                regions.Exon(start=16290, stop=16448),
                regions.Exon(start=16541, stop=16738),
                regions.Exon(start=16916, stop=17427),
                regions.Exon(start=17604, stop=17750),
                regions.Exon(start=17957, stop=18797),
            ),
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
        was_mapped=False,
        identity=None,
        metadata={
            'rna_type': 'lncRNA',
            'providing_databases': ['LNCipedia'],
            'databases': ['LNCipedia', 'NONCODE'],
        }
    ))


def test_coordinates_do_not_exceed_bounds():
    found = fetch_one('URS00008BF974_9606', 'GRCh38')
    assert attr.asdict(found) == attr.asdict(data.Region(
        region_id='URS00008BF974_9606.0',
        rna_id='URS00008BF974_9606',
        region=regions.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='X',
            strand=1,
            exons=(
                regions.Exon(start=114044718, stop=114045067),
                regions.Exon(start=114044725, stop=114044793),
            ),
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
        was_mapped=False,
        identity=None,
        metadata={
            'rna_type': 'lncRNA',
            'providing_databases': ['LNCipedia'],
            'databases': ['LNCipedia'],
        }
    ))


@pytest.mark.xfail()
def test_does_assign_source_correctly():
    found = fetch_one('URS0000000098_9606', 'GRCh38')
    assert attr.asdict(found) == attr.asdict(data.Region(
        region_id='URS0000000098_9606.0',
        rna_id='URS0000000098_9606',
        region=regions.SequenceRegion(
            assembly_id='GRCh38',
            chromosome='10',
            strand=1,
            exons=(regions.Exon(start=17403508, stop=17403618),),
            coordinate_system=regions.CoordinateSystem.one_based(),
        ),
        was_mapped=True,
        identity=1.0,
        metadata={
            'rna_type': 'Y_RNA',
            'databases': ['ENA'],
        },
    ))
