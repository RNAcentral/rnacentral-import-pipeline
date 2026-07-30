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
