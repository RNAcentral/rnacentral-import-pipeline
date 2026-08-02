# -*- coding: utf-8 -*-

"""
Regression tests for the release functions that bulk-insert into FK-bearing
tables.

Postgres queues one FK-check after-trigger event per inserted row and only
drains that queue at end of statement, so a single unbatched INSERT over a full
import dies with::

    psycopg2.errors.OutOfMemory: out of memory
    DETAIL:  Failed on request of size 1048576 in memory context "AfterTriggerEvents".

update_rnc_accessions was fixed for this; update_literature_references was not,
and OOMed on rnc_reference_map after ~2.5 hours. These tests pin the batching
structure so neither can silently regress to one giant statement.
"""

from pathlib import Path

import pytest

FUNCTIONS = Path(__file__).resolve().parents[3] / "database_functions" / "rnc_update"

BATCHED = {
    "update_rnc_accessions.sql": "rnc_accessions",
    "update_literature_references.sql": "rnc_reference_map",
}


@pytest.mark.parametrize("filename,target", sorted(BATCHED.items()))
def test_bulk_insert_is_range_batched(filename, target):
    sql = (FUNCTIONS / filename).read_text()

    # A declared batch size, sliced over a dense rn key in a loop.
    assert "v_batch_size" in sql
    assert "WHILE lo <= v_total LOOP" in sql
    assert "row_number() OVER" in sql
    assert "rn >= %s and rn < %s" in sql.lower()

    # The insert into the FK-bearing table reads from the staging table, not
    # straight from the load table -- the latter is the unbatched shape.
    insert = sql.lower().split(f"insert into rnacen.{target}", 1)
    assert len(insert) == 2, f"no batched insert into rnacen.{target}"
    assert "from load_rnc_references t3" not in insert[1]


def test_reference_map_insert_skips_pairs_that_already_exist():
    """The anti-join is what makes a re-import cheap; ON CONFLICT alone probed
    the index for every pair, almost all of which already existed."""
    sql = (FUNCTIONS / "update_literature_references.sql").read_text().lower()

    assert "left join rnacen.rnc_reference_map m" in sql
    assert "where m.accession is null" in sql
    # Kept as a safety net: the anti-join is not atomic with the insert.
    assert "on conflict (accession, reference_id)" in sql
