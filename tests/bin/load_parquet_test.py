# -*- coding: utf-8 -*-

"""
Tests for the parts of ``bin/load-parquet`` that need no Postgres.

The column-discovery step is the one worth pinning: it used to use
``parquet_schema()``, which returns one row per column per *file*. Every real
caller passes a multi-file glob ('raw*.parquet' in workflows/load-data.nf), so
the generated column list repeated each column once per file and the INSERT
died with "column specified more than once". The bug was invisible while
writer_format defaulted to csv.

Loading the target itself is covered by the shakedown, not here -- it needs a
live staging table.
"""

import importlib.util
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

BIN = Path(__file__).resolve().parents[2] / "bin" / "load-parquet"


def _load_script():
    """Import bin/load-parquet, which has no .py suffix, as a module."""
    spec = importlib.util.spec_from_loader(
        "load_parquet", importlib.machinery.SourceFileLoader("load_parquet", str(BIN))
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


load_parquet = _load_script()


@pytest.fixture
def con():
    connection = duckdb.connect()
    yield connection
    connection.close()


def _write(path, table):
    pq.write_table(table, path)


def test_single_file_columns(con, tmp_path):
    _write(tmp_path / "raw1.parquet", pa.table({"a": [1], "b": ["x"]}))
    columns = load_parquet.parquet_columns(con, str(tmp_path / "raw*.parquet"))
    assert columns == ["a", "b"]


def test_multi_file_glob_does_not_duplicate_columns(con, tmp_path):
    """The regression: two files must still yield one column list."""
    _write(tmp_path / "raw1.parquet", pa.table({"a": [1], "b": ["x"]}))
    _write(tmp_path / "raw2.parquet", pa.table({"a": [2], "b": ["y"]}))

    columns = load_parquet.parquet_columns(con, str(tmp_path / "raw*.parquet"))

    assert columns == ["a", "b"]
    assert len(columns) == len(set(columns))


def test_many_file_glob_does_not_duplicate_columns(con, tmp_path):
    for i in range(5):
        _write(tmp_path / f"raw{i}.parquet", pa.table({"a": [i], "b": ["x"]}))

    columns = load_parquet.parquet_columns(con, str(tmp_path / "raw*.parquet"))

    assert columns == ["a", "b"]


def test_columns_are_unioned_by_name(con, tmp_path):
    """
    Chunks from different parsers can disagree on columns; the column list must
    be the union, in a stable order, still without repeats. This matches the
    ``union_by_name = true`` used by the INSERT ... SELECT.
    """
    _write(tmp_path / "raw1.parquet", pa.table({"a": [1], "b": ["x"]}))
    _write(tmp_path / "raw2.parquet", pa.table({"a": [2], "c": [3.5]}))

    columns = load_parquet.parquet_columns(con, str(tmp_path / "raw*.parquet"))

    assert sorted(columns) == ["a", "b", "c"]
    assert len(columns) == len(set(columns))


def test_generated_column_list_is_valid_sql(con, tmp_path):
    """
    End-to-end within DuckDB: the column list must actually work in the
    INSERT-shaped statement the script builds. A duplicated list fails here
    the same way Postgres would.
    """
    _write(tmp_path / "raw1.parquet", pa.table({"a": [1], "b": ["x"]}))
    _write(tmp_path / "raw2.parquet", pa.table({"a": [2], "b": ["y"]}))
    glob = str(tmp_path / "raw*.parquet")

    columns = load_parquet.parquet_columns(con, glob)
    col_list = ", ".join(f'"{col}"' for col in columns)

    con.execute("CREATE TABLE target (a BIGINT, b VARCHAR)")
    con.execute(
        f"INSERT INTO target ({col_list}) "
        f"SELECT {col_list} FROM read_parquet(?, union_by_name = true)",
        [glob],
    )

    assert con.execute("SELECT count(*) FROM target").fetchone()[0] == 2


def test_reserved_word_columns_are_quoted(con, tmp_path):
    """Column names that are SQL keywords must survive quoting."""
    _write(tmp_path / "raw1.parquet", pa.table({"database": ["ENA"], "order": [1]}))
    glob = str(tmp_path / "raw*.parquet")

    columns = load_parquet.parquet_columns(con, glob)
    col_list = ", ".join(f'"{col}"' for col in columns)

    con.execute('CREATE TABLE target ("database" VARCHAR, "order" BIGINT)')
    con.execute(
        f"INSERT INTO target ({col_list}) "
        f"SELECT {col_list} FROM read_parquet(?, union_by_name = true)",
        [glob],
    )

    assert con.execute("SELECT count(*) FROM target").fetchone()[0] == 1


def test_resource_settings_bound_duckdb_to_the_task():
    """
    Nothing capped DuckDB, so it took hardware_concurrency() threads and ~80% of
    system RAM regardless of the `cpus 4`/`memory 9.GB` the task declared --
    oversubscribing its cores and free to exceed its cgroup.
    """
    settings = load_parquet._resource_settings(4, "8GB")

    assert "SET threads = 4" in settings
    assert "SET memory_limit = '8GB'" in settings
    assert "SET preserve_insertion_order = false" in settings


def test_resource_settings_omits_unset_limits():
    """Bare `bin/load-parquet` runs outside Nextflow keep DuckDB's defaults."""
    settings = load_parquet._resource_settings(None, None)

    assert settings == ["SET preserve_insertion_order = false"]


def test_resource_settings_are_valid_duckdb(con):
    """Each statement has to actually parse and apply."""
    for statement in load_parquet._resource_settings(2, "1GB"):
        con.execute(statement)

    assert con.execute("SELECT current_setting('threads')").fetchone()[0] == 2


class _FakeCursor:
    def __init__(self, recorder, result):
        self.recorder = recorder
        self.result = result

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def execute(self, sql):
        self.recorder.append(sql)

    def fetchone(self):
        return (self.result,)


class _FakeConnection:
    def __init__(self, recorder, result):
        self.recorder = recorder
        self.result = result
        self.closed = False

    def cursor(self):
        return _FakeCursor(self.recorder, self.result)

    def close(self):
        self.closed = True


def test_target_row_count_goes_through_postgres_not_duckdb(monkeypatch):
    """
    The regression: this count used to run over the DuckDB postgres attachment,
    which has no aggregate pushdown and answers it by COPYing the entire table
    back. On a freshly loaded ENA staging table that killed the server with
    'failed to initialize hash table "CompactCheckpointerRequestQueue"' -- after
    the INSERT had committed.
    """
    executed = []
    conn = _FakeConnection(executed, 42)
    monkeypatch.setattr(load_parquet.psycopg2, "connect", lambda url: conn)

    total = load_parquet.target_row_count("postgresql://x", '"rnacen"."load_x"')

    assert total == 42
    assert executed == ['SELECT count(*) FROM "rnacen"."load_x"']
    assert "pg." not in executed[0]
    assert conn.closed


def test_target_row_count_survives_a_database_failure(monkeypatch, capsys):
    """
    The count is only used for a log line, so a hiccup here must not throw away
    a load that already committed.
    """

    def boom(_url):
        raise RuntimeError("server closed the connection unexpectedly")

    monkeypatch.setattr(load_parquet.psycopg2, "connect", boom)

    assert load_parquet.target_row_count("postgresql://x", '"rnacen"."load_x"') is None
    assert "could not count" in capsys.readouterr().err


def test_load_table_mapping_covers_every_writer_output():
    """
    Every key of ENTRY_WRITER_SCHEMAS needs an entry in ENTRY_WRITER_LOAD_TABLES
    -- including an explicit None for outputs with no load target. A missing key
    would be treated as a literal table name and try to load into a table named
    e.g. 'terms'.
    """
    from rnacentral_pipeline.schemas import (
        ENTRY_WRITER_LOAD_TABLES,
        ENTRY_WRITER_SCHEMAS,
    )

    missing = set(ENTRY_WRITER_SCHEMAS) - set(ENTRY_WRITER_LOAD_TABLES)
    assert not missing, f"no load target declared for: {sorted(missing)}"


def test_truncate_issues_a_real_truncate(monkeypatch):
    """
    --truncate must go through psycopg2, not DuckDB's postgres extension.

    The extension rewrites TRUNCATE into DELETE, which leaves the dead tuples
    (and the index bloat) behind -- load_traveler_attempted grew a 1.5 GB pkey
    over 86k rows that way. This pins the statement actually sent to Postgres.
    """
    executed = []

    class FakeCursor:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def execute(self, sql):
            executed.append(sql)

    class FakeConn:
        closed = False

        def set_session(self, autocommit):
            executed.append(f"autocommit={autocommit}")

        def cursor(self):
            return FakeCursor()

        def close(self):
            FakeConn.closed = True

    monkeypatch.setattr(load_parquet.psycopg2, "connect", lambda url: FakeConn())

    load_parquet._truncate("postgresql://ignored", '"rnacen"."load_x"')

    assert executed == ["autocommit=True", 'TRUNCATE "rnacen"."load_x"']
    assert FakeConn.closed, "connection must be closed before the DuckDB INSERT"


def test_count_file_records_rows_loaded(tmp_path):
    """
    load-data.nf reads this file to decide whether the name goes on to select
    pre/post-release SQL, so the value has to be parseable as an int.
    """
    path = tmp_path / "rows.count"
    load_parquet.write_count(path, 5375353)
    assert int(path.read_text().strip()) == 5375353


def test_count_file_records_an_empty_load(tmp_path):
    path = tmp_path / "rows.count"
    load_parquet.write_count(path, 0)
    assert int(path.read_text().strip()) == 0


def test_no_count_file_is_written_when_not_asked(tmp_path):
    """Bare CLI use passes no --count-file and must not litter the cwd."""
    load_parquet.write_count(None, 3)
    assert list(tmp_path.iterdir()) == []
