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

import os
import json

import psycopg2
import pytest

from rnacentral_pipeline.rnacentral.release import database_stats as stats


@pytest.fixture(scope="module")
def connection():
    db_url = os.environ["PGDATABASE"]
    yield psycopg2.connect(db_url)


@pytest.mark.parametrize(
    "db_id",
    [
        4,
    ],
)
def test_gets_correct_length_counts(connection, db_id):
    data = stats.length_counts(connection, db_id)
    with open(f"data/release/lengths-{db_id}.json") as raw:
        expected = json.loads(raw.readline())
    assert json.loads(data) == expected


@pytest.mark.parametrize(
    "db_id",
    [
        # 4,
        28,
        pytest.param(24, marks=pytest.mark.xfail()),
    ],
)
def test_gets_correct_lineage_info(connection, db_id):
    data = stats.lineage(connection, db_id)
    with open(f"data/release/lineage-{db_id}.json") as raw:
        expected = json.loads(raw.readline())
    assert json.loads(data) == expected


@pytest.mark.parametrize(
    "name,expected",
    [
        ("BOB", False),
        ("ENA", True),
    ],
)
def test_knows_it_has_stats(connection, name, expected):
    assert stats.has_stats_for(connection, name) == expected
