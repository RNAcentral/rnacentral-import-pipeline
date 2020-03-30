# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import csv
import io

import pytest

from rnacentral_pipeline.rnacentral import attempted


@pytest.fixture()
def genome_mapping():
    with open('data/genome-mapping/raw.json', 'r') as raw:
        temp = io.StringIO()
        attempted.genome_mapping(raw, 'HanXRQr1.0', temp)
        temp.flush()
        temp.seek(0)
        data = csv.reader(temp)
        yield list(data)

def test_can_get_rfam_version():
    with open('data/attempted/rfam/readme', 'r') as raw:
        assert attempted.parse_rfam_version(raw) == '14.1'


def test_can_produce_all_genome_mapping(genome_mapping):
    assert len(genome_mapping) == 10


def test_can_produce_valid_genome_mapping_entries(genome_mapping):
    assert genome_mapping[0] == ['URS00004AE028_4232', 'HanXRQr1.0']
