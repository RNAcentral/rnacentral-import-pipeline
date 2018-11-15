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

from rnacentral_pipeline.databases.dfam import results


@pytest.fixture(scope='module')
def test_data():
    with open('data/dfam/dfam.tsv', 'r') as raw:
        return list(results.parse(raw))


def test_can_parse_all_results(test_data):
    assert len(test_data) == 8


def test_can_parse_results_correctly(test_data):
    assert attr.asdict(test_data[0]) == attr.asdict(results.Result(
        model_name='MLT1A',
        model_accession='DF0001126.4',
        upi='URS00007CC20E',
        bits=205.7,
        e_value=3.6e-60,
        bias=4.1,
        model_start=1,
        model_end=374,
        strand=1,
        alignment_start=726,
        alignment_end=1062,
    ))
