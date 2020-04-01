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


def write_and_parse(writer, filename, *args):
    with open(filename, 'r') as raw:
        temp = io.StringIO()
        writer(raw, *args, temp)
        temp.flush()
        temp.seek(0)
        data = csv.reader(temp)
        return list(data)


def test_can_get_rfam_version():
    with open('data/attempted/rfam/readme', 'r') as raw:
        assert attempted.parse_rfam_version(raw) == '14.1'


@pytest.mark.parametrize('parser,filename,rest,count', [
    [attempted.genome_mapping, 'data/genome-mapping/raw.json', ['HanXRQr1.0'], 10],
    [attempted.traveler, 'data/traveler/rfam/RF00162/URS0000A7635A.fasta', [], 1],
])
def test_can_produce_all_expected_data(parser, filename, rest, count):
    val = write_and_parse(parser, filename, *rest)
    assert len(val) == count


def test_can_produce_valid_genome_mapping_entries():
    val = write_and_parse(attempted.genome_mapping, 'data/genome-mapping/raw.json', 'HanXRQr1.0')
    assert val[0] == ['URS00004AE028_4232', 'HanXRQr1.0']


def test_can_produce_correct_data_for_traveler():
    val = write_and_parse(attempted.traveler, 'data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.fasta')
    assert val[0] == ['URS00000F9D45_9606']
