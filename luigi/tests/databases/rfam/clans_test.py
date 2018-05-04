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

from rnacentral.mysql import MysqlWrapper

from tasks.config import rfam

from databases.rfam import clans


@pytest.fixture
def stored():
    with open('data/rfam/clans.tsv', 'r') as raw:
        yield raw


@pytest.fixture
def computed():
    with tempfile.NamedTemporaryFile() as tmp:
        mysql = MysqlWrapper(rfam())
        mysql.write_query(tmp, clans.QUERY)
        tmp.flush()
        with open(tmp.name, 'r') as raw:
            yield raw



@pytest.mark.skip()
def test_can_parse_known_data(stored):
    pass


@pytest.mark.skip()
def test_can_parse_queried_data(computed):
    pass


@pytest.mark.skip()
def test_can_prdouce_correct_values(stored):
    pass
