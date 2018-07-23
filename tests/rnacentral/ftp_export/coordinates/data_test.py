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

import os

import pytest

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import data

from tests.helpers import run_with_replacements


def fetch_raw(rna_id, assembly):
    path = os.path.join('files', 'ftp-export', 'genome_coordinates',
                        'query.sql')
    _, taxid = rna_id.split('_')
    return run_with_replacements(
        path,
        (':assembly_id', "'GRCh38'"),
        (':taxid', taxid),
        ('WHERE\n',
         '''where
         pre.id = '%s'
         and''' % rna_id
         ), take_all=True)


def fetch_one(rna_id, assembly):
    return next(data.parse(fetch_raw(rna_id, assembly)))


@pytest.mark.parametrize('rna_id,assembly,count', [
    ('URS0000A78C33_9606', 'GRCh38', 4),
])
def test_can_fetch_all_coordinates_for_upi_taxid(rna_id, assembly, count):
    assert len(fetch_raw(rna_id, assembly)) == count


def test_can_build_correct_data_for_known():
    found = fetch_one('URS0000A78C33_9606', 'GRCh38')

    assert found.regions[0].start == 5126
    assert found.regions[0].stop == 50017
    assert found == data.LocatedSequence(
        rna_id='URS0000A78C33_9606',
        rna_type='lncRNA',
        databases="Ensembl,GENCODE",
        regions=[
            data.Region(
                region_id='URS0000A78C33_9606.1',
                chromosome='CHR_HSCHR5_1_CTG1',
                strand=1,
                endpoints=(
                    data.Endpoint(start=5126, stop=5218),
                    data.Endpoint(start=49916, stop=50017),
                ),
            ),
        ],
    )


def test_can_build_correct_data_for_both_mapped_and_known():
    found = fetch_one('URS00008C1902_9606', 'GRCh38')

    assert found.regions[0].start == 11869
    assert found.regions[0].stop == 14409
    import attr
    from pprint import pprint
    pprint(attr.asdict(found))
    assert found == data.LocatedSequence(
        rna_id='URS00008C1902_9606',
        rna_type='lncRNA',
        databases="LNCipedia,NONCODE",
        regions=[
            data.Region(
                region_id='URS00008C1902_9606.1',
                chromosome='1',
                strand=1,
                endpoints=(
                    data.Endpoint(start=11869, stop=12227),
                    data.Endpoint(start=12613, stop=12721),
                    data.Endpoint(start=13221, stop=14409),
                ),
            ),
        ],
    )
