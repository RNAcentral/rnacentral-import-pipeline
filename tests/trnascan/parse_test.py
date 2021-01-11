# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from pathlib import Path

import pytest

from rnacentral_pipeline.databases.data import TRnaScanResults
from rnacentral_pipeline.trnascan import parser


@pytest.mark.parametrize(
    "filename,count",
    [
        ("data/r2dt/gtrnadb/E-shuffled.part_2464.txt", 223),
        ("data/r2dt/gtrnadb/B-shuffled.part_2464.txt", 216),
        ("data/r2dt/gtrnadb/A-shuffled.part_2464.txt", 166),
    ],
)
def test_can_parse_a_file(filename, count):
    path = Path(filename)
    assert len(list(parser.parse(path))) == count


def test_produces_reasonable_data():
    path = Path("data/r2dt/gtrnadb/E-shuffled.part_2464.txt")
    assert next(parser.parse(path)) == TRnaScanResults(
        sequence_id='URS0000C7FBE7',
        hit_index=1,
        sequence_start=1,
        sequence_stop=73,
        trna_type='Val',
        anticodon='TAC',
        intron_start=None,
        intron_stop=None,
        score=40.6,
        note='pseudo',
    )
