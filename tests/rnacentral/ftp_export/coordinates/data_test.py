# -*- coding: utf-8 -*-

# pylint: disable=no-member

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

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data

from .helpers import fetch_coord, fetch_all


def fetch_one(rna_id, assembly):
    coords = list(fetch_coord(rna_id, assembly))
    assert len(coords) == 1
    return coords[0]


@pytest.mark.parametrize('rna_id,assembly,count', [
    # ('URS0000A78C33_9606', 'GRCh38', 4),  # Not imported?
    # ('URS00009BF201_9606', 'GRCh38', 9),  # What to do with mapped and given?
    ('URS00008B37EC_9606', 'GRCh38', 1),
    ('URS00008C1914_9606', 'GRCh38', 1),
    ('URS00006683B1_281687', 'GRCh38', 0),
])
def test_can_fetch_all_coordinates_for_upi_taxid(rna_id, assembly, count):
    assert len(list(fetch_coord(rna_id, assembly))) == count


@pytest.mark.skip
@pytest.mark.parametrize('assembly,count', [  # pylint: disable=no-member
    ('GRCh38', 383001),
    ('BDGP6', 40602),
    ('WBcel235', 29910),
    ('TAIR10', 277777),
    ('GRCm38', 978299),
    ('ASM294v2', 3310),
    ('SL2.50', 4677),
])
def test_can_find_all_required_coordinates(assembly, count):
    assert len(list(fetch_all(assembly))) == count


@pytest.mark.skip
def test_can_build_for_only_mapped():
    found = fetch_one('URS000012C1C6_9606', 'GRCh38')

    ans = [
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='1',
            strand=1,
            endpoints=(data.Endpoint(start=165777430, stop=165777459),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='3',
            strand=-1,
            endpoints=(data.Endpoint(start=163618677, stop=163618706),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='8',
            strand=-1,
            endpoints=(data.Endpoint(start=52338512, stop=52338541),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='10',
            strand=-1,
            endpoints=(data.Endpoint(start=73941772, stop=73941801),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='12',
            strand=1,
            endpoints=(data.Endpoint(start=123897625, stop=123897654),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='19',
            strand=-1,
            endpoints=(data.Endpoint(start=3161031, stop=3161060),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='22',
            strand=-1,
            endpoints=(data.Endpoint(start=18174677, stop=18174706),),
            source='alignment',
            identity=1.0
        ),
        data.Region(
            rna_id='URS000012C1C6_9606',
            rna_type='piRNA',
            databases=['ENA'],
            chromosome='22',
            strand=1,
            endpoints=(data.Endpoint(start=20362774, stop=20362803),),
            source='alignment',
            identity=1.0
        ),
    ]

    assert attr.asdict(found) == attr.asdict(ans)


@pytest.mark.skip
def test_can_find_correct_for_something_that_can_be_mapped():
    found = fetch_one('URS00009BF201_9606', 'GRCh38')
    assert found == [
        data.Region(
            rna_id='URS00009BF201_9606',
            rna_type='lncRNA',
            databases=['LNCipedia', 'NONCODE'],
            chromosome='16',
            strand=-1,
            endpoints=(
                data.Endpoint(start=14085, stop=14511),
                data.Endpoint(start=14652, stop=14720),
                data.Endpoint(start=15481, stop=15633),
                data.Endpoint(start=16290, stop=16448),
                data.Endpoint(start=16541, stop=16738),
                data.Endpoint(start=16916, stop=17427),
                data.Endpoint(start=17604, stop=17750),
                data.Endpoint(start=17957, stop=18797),
            )
        )
    ]
