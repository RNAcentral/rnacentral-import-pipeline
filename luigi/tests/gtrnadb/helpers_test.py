# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

from databases.gtrnadb import helpers


@pytest.fixture
def data():
    with open('data/gtrnadb.json', 'rb') as raw:
        return json.load(raw)


def test_url(data):
    assert helpers.url(data[0]) == "http://gtrnadb.ucsc.edu/genomes/bacteria/Acar_mari_MBIC11017/genes/tRNA-Ala-CGC-1-1.html"


def test_anticodon(data):
    assert helpers.anticodon(data[0]) == 'CGC'


def test_note_data(data):
    assert helpers.note_data(data[0]) == {
        "anticodon": 'CGC',
        "anticodon_positions": [
            {
                "relative_start": 34,
                "relative_stop": 36
            }
        ],
        "isotype": "Ala",
        "score": 72.7,
        "url": "http://gtrnadb.ucsc.edu/genomes/bacteria/Acar_mari_MBIC11017/genes/tRNA-Ala-CGC-1-1.html"
    }


def test_common_name(data):
    assert helpers.common_name(data[0]) is None


def test_lineage(data):
    assert helpers.lineage(data[0]) == 'Bacteria; Cyanobacteria; Synechococcales; Acaryochloridaceae; Acaryochloris; Acaryochloris marina MBIC11017'


def test_species(data):
    assert helpers.species(data[0]) == "Acaryochloris marina MBIC11017"


def test_product(data):
    assert helpers.product(data[0]) == 'tRNA-Ala (CGC)'


def test_as_dotbracket(data):
    ans = "(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))."
    assert helpers.dot_bracket(data[0]) == ans


def test_simple_description(data):
    assert helpers.description(data[0]) == "Acaryochloris marina MBIC11017 tRNA-Ala (CGC)"


def test_adds_label_about_introns(data):
    entry = next(d for d in data if d['gene'] == 'tRNA-Ser-CGA-2-1')
    name = 'Acaryochloris marina MBIC11017 tRNA-Ser (CGA) containing a group I intron'
    assert helpers.description(entry) == name


def test_as_dotbracket_detects_weird_strings():
    data = {'secondary_structure': '>>>...A<<<'}
    with pytest.raises(helpers.InvalidDotBracket):
        helpers.dot_bracket(data)


def test_primary_id_is_always_unique(data):
    seen = set()
    for entry in data:
        for location in entry['genome_locations']:
            pid = helpers.primary_id(entry, location)
            assert pid not in seen
            seen.add(pid)


def test_builds_primary_id(data):
    pids = []
    entry = data[0]
    for location in entry['genome_locations']:
        pid = helpers.primary_id(entry, location)
        pids.append(pid)
    assert pids == [
        "tRNA-Ala-CGC-1-1:CP000828.1:603738-603810"
    ]


def test_chromosome(data):
    assert helpers.chromosome(data[0]['genome_locations'][0]) == 'chr'
