# -*- coding: utf-8 -*-

"""
Tests for the database function deployer.

No Postgres: ``apply()`` is driven through a recording fake cursor, which is
enough to pin the ordering property that mattered -- the tracking table has to
be created before anything reads it. Nothing in the repo creates
``rnacen.applied_functions`` (no DDL under database_functions/, no migration),
so without the bootstrap the first run against any database died on the very
first SELECT.

``discover()`` and ``pending()`` are pure and run against the real
database_functions/ tree.
"""

import hashlib
from pathlib import Path

import pytest

from rnacentral_pipeline.rnacentral.release import functions

TRACKING = functions.TRACKING_TABLE


class FakeCursor:
    """Records executed SQL. Returns no applied rows unless told otherwise."""

    def __init__(self, applied=()):
        self.executed = []
        self._applied = list(applied)

    def execute(self, sql, params=None):
        self.executed.append((" ".join(sql.split()), params))

    def fetchall(self):
        return self._applied

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class FakeConnection:
    def __init__(self, cursor):
        self._cursor = cursor
        self.committed = False
        self.rolled_back = False

    def cursor(self):
        return self._cursor

    def commit(self):
        self.committed = True

    def rollback(self):
        self.rolled_back = True

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


@pytest.fixture
def fake_db(monkeypatch):
    cursor = FakeCursor()
    conn = FakeConnection(cursor)
    monkeypatch.setattr(functions.psycopg2, "connect", lambda url: conn)
    return conn, cursor


def _sql_at(cursor, index):
    return cursor.executed[index][0]


# ---------------------------------------------------------------------------
# discover()


def test_discover_finds_real_functions():
    found = functions.discover()
    assert found
    assert all(f.path.suffix == ".sql" for f in found)
    assert all(f.sql for f in found)


def test_discover_is_sorted_by_schema_then_name():
    found = functions.discover()
    keys = [f.key for f in found]
    assert keys == sorted(keys)


def test_discover_excludes_requested_schemas():
    everything = functions.discover()
    assert any(f.schema == "rnc_test" for f in everything), (
        "expected rnc_test fixtures in database_functions/; "
        "the exclusion test below is meaningless without them"
    )

    deployable = functions.discover(exclude_schemas={"rnc_test"})
    assert not any(f.schema == "rnc_test" for f in deployable)
    assert len(deployable) < len(everything)


def test_deploy_excludes_rnc_test_by_default():
    """rnc_test holds setup/teardown fixtures that must never hit production."""
    assert "rnc_test" in functions.DEPLOY_EXCLUDE_SCHEMAS


def test_discover_raises_on_empty_directory(tmp_path):
    with pytest.raises(ValueError):
        functions.discover(base=tmp_path)


def test_sha256_tracks_file_contents(tmp_path):
    fn = functions.Function(
        schema="rnacen", name="thing", path=tmp_path / "thing.sql", sql="SELECT 1"
    )
    assert fn.sha256 == hashlib.sha256(b"SELECT 1").hexdigest()


# ---------------------------------------------------------------------------
# pending()


def _function(schema, name, sql):
    return functions.Function(
        schema=schema, name=name, path=Path(f"{schema}/{name}.sql"), sql=sql
    )


def test_pending_returns_everything_when_nothing_applied():
    fns = [_function("rnacen", "a", "SELECT 1")]
    assert functions.pending(FakeCursor(), fns) == fns


def test_pending_skips_unchanged_functions():
    fn = _function("rnacen", "a", "SELECT 1")
    cursor = FakeCursor(applied=[("rnacen", "a", fn.sha256)])
    assert functions.pending(cursor, [fn]) == []


def test_pending_detects_changed_sha():
    fn = _function("rnacen", "a", "SELECT 1")
    cursor = FakeCursor(applied=[("rnacen", "a", "stale-sha")])
    assert functions.pending(cursor, [fn]) == [fn]


# ---------------------------------------------------------------------------
# apply() -- the bootstrap


def test_apply_reads_the_tracking_table_it_never_creates(fake_db):
    """
    apply() assumes rnacen.applied_functions already exists -- nothing in the
    repo creates it. Documented in
    docs/database-functions-vault/applied-functions-tracking-table.md; pinned
    here so the assumption is visible rather than implicit.
    """
    _, cursor = fake_db
    functions.apply("postgresql://ignored")

    assert not any("CREATE TABLE" in sql for sql, _ in cursor.executed)
    assert any(
        sql.startswith(f"SELECT") and TRACKING in sql for sql, _ in cursor.executed
    )


def test_apply_records_each_applied_function(fake_db):
    conn, cursor = fake_db
    applied = functions.apply("postgresql://ignored")

    assert applied
    inserts = [
        params
        for sql, params in cursor.executed
        if sql.startswith(f"INSERT INTO {TRACKING}")
    ]
    assert len(inserts) == len(applied)
    assert {(p[0], p[1]) for p in inserts} == {f.key for f in applied}
    assert conn.committed


def test_apply_commits_only_once_at_the_end(fake_db):
    """All-or-nothing: one transaction covering every function."""
    conn, _ = fake_db
    functions.apply("postgresql://ignored")
    assert conn.committed
    assert not conn.rolled_back


def test_dry_run_creates_nothing_and_rolls_back(fake_db):
    conn, cursor = fake_db
    would_apply = functions.apply("postgresql://ignored", dry_run=True)

    assert would_apply
    assert conn.rolled_back
    assert not conn.committed
    assert not any(
        sql.startswith(f"INSERT INTO {TRACKING}") for sql, _ in cursor.executed
    )


def test_apply_is_a_noop_when_everything_is_current(monkeypatch):
    fns = functions.discover(exclude_schemas=functions.DEPLOY_EXCLUDE_SCHEMAS)
    cursor = FakeCursor(applied=[(f.schema, f.name, f.sha256) for f in fns])
    conn = FakeConnection(cursor)
    monkeypatch.setattr(functions.psycopg2, "connect", lambda url: conn)

    assert functions.apply("postgresql://ignored") == []
    assert not any(
        sql.startswith(f"INSERT INTO {TRACKING}") for sql, _ in cursor.executed
    )


def test_apply_never_touches_excluded_schemas(fake_db):
    applied = functions.apply("postgresql://ignored")
    assert not any(f.schema in functions.DEPLOY_EXCLUDE_SCHEMAS for f in applied)
