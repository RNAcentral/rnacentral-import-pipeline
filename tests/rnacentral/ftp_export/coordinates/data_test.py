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

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data

from .helpers import fetch_raw, fetch_coord


def fetch_one(rna_id, assembly):
    return next(fetch_coord(rna_id, assembly))


@pytest.mark.parametrize('rna_id,assembly,count', [
    ('URS0000A78C33_9606', 'GRCh38', 4),
])
def test_can_fetch_all_coordinates_for_upi_taxid(rna_id, assembly, count):
    assert len(fetch_raw(rna_id, assembly)) == count


def test_can_build_correct_data_for_known():
    found = fetch_one('URS0000A78C33_9606', 'GRCh38')

    assert found.regions[0].start == 5126
    assert found.regions[0].stop == 50017
    assert found.regions[0].region_id == 'URS0000A78C33_9606@CHR_HSCHR5_1_CTG1/5126-50017:+'
    assert attr.asdict(found) == attr.asdict(data.LocatedSequence(
        rna_id='URS0000A78C33_9606',
        rna_type='lncRNA',
        databases="Ensembl,GENCODE",
        regions=[
            data.Region(
                rna_id='URS0000A78C33_9606',
                chromosome='CHR_HSCHR5_1_CTG1',
                strand=1,
                endpoints=(
                    data.Endpoint(start=5126, stop=5218),
                    data.Endpoint(start=49916, stop=50017),
                ),
                source='expert-database',
            ),
        ],
    ))


def test_can_build_correct_data_for_both_mapped_and_known():
    found = fetch_one('URS00008C1902_9606', 'GRCh38')

    ans = data.LocatedSequence(
        rna_id='URS00008C1902_9606',
        rna_type='lncRNA',
        databases="LNCipedia,NONCODE",
        regions=[
            data.Region(
                rna_id='URS00008C1902_9606',
                chromosome='1',
                strand=1,
                endpoints=(
                    data.Endpoint(start=11869, stop=12227),
                    data.Endpoint(start=12613, stop=12721),
                    data.Endpoint(start=13221, stop=14409),
                ),
                source='expert-database',
            ),
        ],
    )

    assert found.regions[0].start == 11869
    assert found.regions[0].stop == 14409
    assert found.regions[0].region_id == 'URS00008C1902_9606@1/11869-14409:+'
    assert attr.asdict(found) == attr.asdict(ans)


def test_can_build_for_only_mapped():
    found = fetch_one('URS000012C1C6_9606', 'GRCh38')

    ans = data.LocatedSequence(
        rna_id='URS000012C1C6_9606',
        rna_type='piRNA',
        databases='ENA',
        regions=[
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='1',
                strand=1,
                endpoints=(data.Endpoint(start=165777430, stop=165777459),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='3',
                strand=-1,
                endpoints=(data.Endpoint(start=163618677, stop=163618706),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='8',
                strand=-1,
                endpoints=(data.Endpoint(start=52338512, stop=52338541),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='10',
                strand=-1,
                endpoints=(data.Endpoint(start=73941772, stop=73941801),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='12',
                strand=1,
                endpoints=(data.Endpoint(start=123897625, stop=123897654),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='19',
                strand=-1,
                endpoints=(data.Endpoint(start=3161031, stop=3161060),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='22',
                strand=-1,
                endpoints=(data.Endpoint(start=18174677, stop=18174706),),
                source='alignment',
                identity=1.0
            ),
            data.Region(
                rna_id='URS000012C1C6_9606',
                chromosome='22',
                strand=1,
                endpoints=(data.Endpoint(start=20362774, stop=20362803),),
                source='alignment',
                identity=1.0
            ),
        ]
    )

    assert attr.asdict(found) == attr.asdict(ans)
