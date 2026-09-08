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

import json
import os

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
    "db_id,count",
    [
        (28, 11124),
    ],
)
def test_counts_sequences_correctly(connection, db_id, count):
    found = stats.count_sequences(connection, db_id)
    assert found == count


@pytest.mark.parametrize(
    "db_id,count",
    [
        (28, 1),
    ],
)
def test_counts_organisms_correctly(connection, db_id, count):
    found = stats.count_organisms(connection, db_id)
    assert found == count


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


@pytest.fixture
def stale_fixture(connection):
    """
    TEMP tables shadow the real rnc_database/release_stats for this session, so
    STALE_DATABASES_QUERY can be exercised on known data and rolled back.
    """
    with connection.cursor() as cur:
        cur.execute(
            """
            CREATE TEMP TABLE rnc_database (
                id bigint, descr text, alive char,
                last_import_date timestamp
            ) ON COMMIT DROP;
            CREATE TEMP TABLE release_stats (
                dbid bigint, end_time timestamp
            ) ON COMMIT DROP;

            INSERT INTO rnc_database VALUES
                (1, 'LOADED',     'Y', '2026-01-01'),  -- loaded since last stats
                (2, 'UNCHANGED',  'Y', '2026-01-01'),  -- stats already current
                (3, 'NEVER',      'Y', '2026-01-01'),  -- no load ever recorded
                (4, 'NEW',        'Y', NULL),          -- no stats yet
                (5, 'DEAD',       'N', '2026-01-01');  -- not alive
            INSERT INTO release_stats VALUES
                (1, '2026-02-01'),
                (2, '2025-12-01'),
                (5, '2026-02-01');
            """
        )
        yield cur
    connection.rollback()


def selected(cur, force):
    cur.execute(stats.STALE_DATABASES_QUERY, {"force": force})
    return {descr for (_, descr) in cur}


def test_only_selects_databases_loaded_since_last_stats(stale_fixture):
    assert selected(stale_fixture, False) == {"LOADED", "NEW"}


def test_force_selects_every_live_database(stale_fixture):
    assert selected(stale_fixture, True) == {"LOADED", "UNCHANGED", "NEVER", "NEW"}
