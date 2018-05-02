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

import tempfile

import pytest

import databases.rfam.cross_references as cr


@pytest.fixture
def data():
    with open('data/rfam/database_link.tsv', 'r') as raw:
        return list(cr.parse(raw))


def test_can_fetch_and_parse_data(data):
    assert len(data) == 7909


def test_correctly_parses_so_data(data):
    assert data[0] == cr.RfamDatabaseLink(
        'RF00014',
        'SO',
        None,
        'SO:0000378',
        'DsrA_RNA',
    )


def test_correctly_parses_other(data):
    assert data[4] == cr.RfamDatabaseLink(
        'RF00016',
        'snoopy',
        None,
        'Mus_musculus300888;',
        None,
    )
