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

from rnacentral_pipeline.databases.crs import parser as crs


def test_can_parse_complete_file():
    with open('data/crs/hg38.tsv', 'r') as raw:
        assert len(list(crs.parse(raw))) == 10007


def test_builds_correct_entries():
    with open('data/crs/hg38.tsv', 'r') as raw:
        data = list(f for f in crs.parse(raw) if f.upi == 'URS00009BF201')

    assert len(data) == 9
    assert attr.asdict(data[0]) == attr.asdict(
        crs.CrsFeature(
            upi='URS00009BF201',
            taxid=9606,
            crs_name='M1412625',
            start=763,
            stop=802,
            genomic_location=crs.GenomicLocation(
                chromosome='1',
                start=18307,
                stop=18346,
            ),
            fdr=24.43,
        )
    )

    assert attr.asdict(data[-1]) == attr.asdict(
        crs.CrsFeature(
            upi='URS00009BF201',
            taxid=9606,
            crs_name='M2543977',
            start=-23,
            stop=72,
            genomic_location=crs.GenomicLocation(
                chromosome='X',
                start=156025217,
                stop=156025312,
            ),
            fdr=7.37,
        )
    )
