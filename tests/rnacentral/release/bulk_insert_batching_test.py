# -*- coding: utf-8 -*-

"""
Regression tests for the release functions that bulk-insert into FK-bearing
tables.

Postgres queues one FK-check after-trigger event per inserted row and only
drains that queue at end of statement, so a single unbatched INSERT over a full
import dies with::

    psycopg2.errors.OutOfMemory: out of memory
    DETAIL:  Failed on request of size 1048576 in memory context "AfterTriggerEvents".

update_rnc_accessions is batched for this reason; the test below pins that
structure so it cannot silently regress to one giant statement.

update_literature_references is deliberately NOT batched: its rnc_reference_map
insert anti-joins against what is already there, so on a re-import it carries
only genuinely new pairs. The anti-join is what keeps it small, and it has its
own test here.
"""

from pathlib import Path

import pytest

FUNCTIONS = Path(__file__).resolve().parents[3] / "database_functions" / "rnc_update"


def test_accessions_bulk_insert_is_range_batched():
    sql = (FUNCTIONS / "update_rnc_accessions.sql").read_text()

    # A declared batch size, sliced over a dense rn key in a loop.
    assert "v_batch_size" in sql
    assert "WHILE lo <= v_total LOOP" in sql
    assert "row_number() OVER" in sql
    assert "rn >= %s and rn < %s" in sql.lower()

    # The insert reads from the deduped staging table, not straight from the
    # load table -- the latter is the unbatched shape.
    insert = sql.lower().split("insert into rnacen.rnc_accessions", 1)
    assert len(insert) == 2, "no batched insert into rnacen.rnc_accessions"
    assert "from load_dedup" in insert[1]
    assert "from rnacen.load_rnc_accessions" not in insert[1]


def test_reference_map_insert_skips_pairs_that_already_exist():
    """The anti-join is what makes a re-import cheap; ON CONFLICT alone probed
    the index for every pair, almost all of which already existed."""
    sql = (FUNCTIONS / "update_literature_references.sql").read_text().lower()

    assert "left join rnacen.rnc_reference_map m" in sql
    assert "where m.accession is null" in sql
    # Kept as a safety net: the anti-join is not atomic with the insert.
    assert "on conflict (accession, reference_id)" in sql
