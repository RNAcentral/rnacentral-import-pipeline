# -*- coding: utf-8 -*-

"""
xref is partitioned by dbid. Postgres only skips partitions when that key is
compared to a constant or a parameter, so joining ``x.dbid = l.in_dbid`` (a column
from the load table) prunes nothing and the plan reads all 59 partitions. That is
what made the delta load of ENA slower than the full loads it replaced.

The same mistake has been fixed before in populate_pel_tables1/2 and
log_release_end_atx; these three steps were missed.
"""

from pathlib import Path

import pytest

FUNCTIONS = Path(__file__).resolve().parents[3] / "database_functions"
INCREMENTAL = FUNCTIONS / "rnc_load_xref_incremental"

# The steps that join or update xref and used to receive no dbid at all.
DBID_SCOPED = [
    "incremental_new_versions",
    "incremental_refresh",
    "incremental_retire_changed",
]

CALLERS = ["load_xref_incremental.sql", "load_xref_delta.sql"]

# How each step is called once it takes the dbid. new_versions and refresh only
# ever write the current release, which is already on every staging row, so they
# take the dbid alone; retire_changed still stamps the previous release.
CALLS = {
    "incremental_new_versions": "incremental_new_versions(p_in_dbid)",
    "incremental_refresh": "incremental_refresh(p_in_dbid)",
    "incremental_retire_changed": "incremental_retire_changed(p_in_dbid, p_previous_release)",
}


@pytest.mark.parametrize("name", DBID_SCOPED)
def test_step_restricts_xref_to_one_partition(name):
    sql = (INCREMENTAL / f"{name}.sql").read_text()

    assert f"{name}(p_in_dbid bigint" in sql, "dbid must be a parameter, not a join key"
    assert "dbid = p_in_dbid" in sql


@pytest.mark.parametrize("name", DBID_SCOPED)
def test_old_one_arg_signature_is_dropped(name):
    """CREATE OR REPLACE cannot change an argument list; it would leave an overload."""
    sql = (INCREMENTAL / f"{name}.sql").read_text()

    assert f"DROP FUNCTION IF EXISTS rnc_load_xref_incremental.{name}(bigint);" in sql


@pytest.mark.parametrize("caller", CALLERS)
def test_callers_pass_the_dbid(caller):
    sql = (INCREMENTAL / caller).read_text()

    for name in DBID_SCOPED:
        assert CALLS[name] in sql, name
