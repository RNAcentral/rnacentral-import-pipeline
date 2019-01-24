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

import pytest

from rnacentral_pipeline.databases.ensembl.metadata import karyotypes as karyo


def karyotype(domain, species):
    raw = karyo.fetch(species, domain)
    return karyo.process(raw)


def test_builds_empty_karyotype_for_missing_data():
    _, found = karyotype('ensemblgenomes', 'glycine_max')
    assert len(found) == 1190
    assert found['1'] == {
        'size': 56831624,
        'bands': [{
            'start': 1,
            'end': 56831624,
        }]
    }


def test_builds_with_known_bands():
    _, found = karyotype('ensembl', 'homo_sapiens')
    assert len(found) == 194
    assert found['MT'] == {
        'size': 16569,
        'bands': [{
            'start': 1,
            'end': 16569,
        }]
    }
    assert found['Y'] == {
        'size': 57227415,
        'bands': [
            {"type": "acen", "id": "p11.1", "start": 10300001, "end": 10400000},
            {"id": "p11.2", "end": 10300000, "start": 600001, "type": "gneg"},
            {"end": 600000, "start": 300001, "id": "p11.31", "type": "gpos50"},
            {"end": 300000, "start": 1, "id": "p11.32", "type": "gneg"},
            {"start": 10400001, "end": 10600000, "id": "q11.1", "type": "acen"},
            {"type": "gneg", "start": 10600001, "end": 12400000, "id": "q11.21"},
            {"id": "q11.221", "end": 17100000, "start": 12400001, "type": "gpos50"},
            {"end": 19600000, "start": 17100001, "id": "q11.222", "type": "gneg"},
            {"type": "gpos50", "id": "q11.223", "end": 23800000, "start": 19600001},
            {"type": "gneg", "id": "q11.23", "end": 26600000, "start": 23800001},
            {"type": "gvar", "id": "q12", "start": 26600001, "end": 57227415}
        ]
    }
