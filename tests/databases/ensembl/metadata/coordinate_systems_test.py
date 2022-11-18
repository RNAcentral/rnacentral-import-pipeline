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

import json

import pytest

from rnacentral_pipeline.databases.ensembl.metadata import coordinate_systems as coord


@pytest.fixture()
def chicken():
    with open("data/ensembl/coordinate-systems/chicken.json") as raw:
        return json.load(raw)


@pytest.fixture()
def human():
    with open("data/ensembl/coordinate-systems/human.json") as raw:
        return json.load(raw)


@pytest.mark.ensembl
@pytest.mark.skip()
def test_can_generate_simple_karyotype_information(human):
    data = coord.top_level_only(human)
    data = {c[0]: c[4] for c in data}
    assert data == {
        "1": 1,
        "2": 2,
        "3": 3,
        "4": 4,
        "5": 5,
        "6": 6,
        "7": 7,
        "8": 8,
        "9": 9,
        "10": 10,
        "11": 11,
        "12": 12,
        "13": 13,
        "14": 14,
        "15": 15,
        "16": 16,
        "17": 17,
        "18": 18,
        "19": 19,
        "20": 20,
        "21": 21,
        "22": 22,
        "X": 23,
        "Y": 24,
        "MT": 25,
    }


def test_can_generate_strange_karyotype_information(chicken):
    data = coord.top_level_only(chicken)
    data = {c[0]: c[4] for c in data}
    assert data == {
        "MT": 0,
        "W": 1,
        "Z": 2,
        "1": 3,
        "2": 4,
        "3": 5,
        "4": 6,
        "5": 7,
        "6": 8,
        "7": 9,
        "8": 10,
        "9": 11,
        "10": 12,
        "11": 13,
        "12": 14,
        "13": 15,
        "14": 16,
        "15": 17,
        "16": 18,
        "17": 19,
        "18": 20,
        "19": 21,
        "20": 22,
        "21": 23,
        "22": 24,
        "23": 25,
        "24": 26,
        "25": 27,
        "26": 28,
        "27": 29,
        "28": 30,
        "30": 31,
        "31": 32,
        "32": 33,
        "33": 34,
    }
