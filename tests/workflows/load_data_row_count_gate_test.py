# -*- coding: utf-8 -*-

"""
A loaded name is what selects a pre/post-release SQL script in load-data.nf, so
whatever decides "this parser produced nothing" also decides which release
scripts run.

Under csv that decision was ``filter { f -> !f.isEmpty() }``: an empty parser
left a 0-byte file and its name never reached the glob. A zero-row parquet file
still carries schema and footer bytes, so the same filter passes it through --
which is how 001__karyotypes.sql, dead since 2019, got selected and killed a
release on a table that does not exist.

The gate is now the row count load-parquet reports. These pin the wiring end to
end, because nothing else fails loudly if it comes apart.
"""

import re
from pathlib import Path

import pyarrow as pa
import pytest

from rnacentral_pipeline.parquet_writers import parquet_writer

ROOT = Path(__file__).resolve().parents[2]
LOAD_DATA = (ROOT / "workflows" / "load-data.nf").read_text()


def test_zero_row_parquet_is_not_empty(tmp_path):
    """The premise: file size cannot tell you a parquet file has no rows."""
    path = tmp_path / "empty.parquet"
    schema = pa.schema([pa.field("a", pa.string())])
    with parquet_writer(path, schema):
        pass

    assert path.stat().st_size > 0


def merge_and_import() -> str:
    body = re.search(r"process merge_and_import \{(.*?)\n\}", LOAD_DATA, re.S)
    assert body, "merge_and_import not found in load-data.nf"
    return body.group(1)


def test_merge_and_import_emits_the_count():
    assert "path('rows.count')" in merge_and_import()


@pytest.mark.parametrize("branch", re.findall(r'"""(.*?)"""', merge_and_import(), re.S))
def test_every_load_branch_writes_the_count(branch):
    """
    Both writer formats reach the same output block, so a branch that skips
    rows.count fails the task on a missing output rather than silently
    ungating the release scripts.
    """
    assert "rows.count" in branch


def test_load_data_gates_imported_names_on_the_count():
    """
    Without this filter the name flows on regardless of the count, and an empty
    parser selects its release scripts again.
    """
    body = re.search(r"workflow load_data \{(.*)", LOAD_DATA, re.S)
    assert body, "workflow load_data not found in load-data.nf"

    gate = re.search(
        r"\| merge_and_import(.*?)set \{ imported_names \}", body.group(1), re.S
    )
    assert gate, "imported_names is no longer built from merge_and_import"
    assert "filter" in gate.group(1)
    assert "toInteger() > 0" in gate.group(1)
