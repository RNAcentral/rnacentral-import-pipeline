# -*- coding: utf-8 -*-

"""
Tests for the SQL that decides whether an R2DT diagram is shown.

The query files are executed verbatim, so these guard the shipped artifacts:

  files/search-export/parts/r2dt.sql     - what the search index gets
  files/precompute/queries/r2dt-hits.sql - what precompute gets
  files/r2dt/should-show/update.ctl      - how classifier output reaches the DB
"""

import io
import json
import re
from pathlib import Path

import pytest

SEARCH_EXPORT = Path("files/search-export/parts/r2dt.sql")
PRECOMPUTE = Path("files/precompute/queries/r2dt-hits.sql")
UPDATE_CTL = Path("files/r2dt/should-show/update.ctl")

pytestmark = [pytest.mark.r2dt, pytest.mark.db]


CTL_SECTIONS = ("BEFORE LOAD DO", "WITH", "AFTER LOAD DO")


def ctl_blocks(section):
    """The $$-quoted SQL statements under a DO section of the pgloader CTL."""
    text = UPDATE_CTL.read_text()
    start = text.index(section) + len(section)
    later = [text.index(s) for s in CTL_SECTIONS if s != section]
    stop = min((p for p in later if p > start), default=len(text))
    blocks = re.findall(r"\$\$(.*?)\$\$", text[start:stop], flags=re.DOTALL)
    assert blocks, f"No SQL found under {section} in {UPDATE_CTL}"
    return [b.strip() for b in blocks]


def copy_query(db, path):
    """Run a `COPY (...) TO STDOUT` file and parse the JSON it emits."""
    out = io.StringIO()
    with db.cursor() as cur:
        cur.copy_expert(path.read_text(), out)
    return [json.loads(line) for line in out.getvalue().splitlines() if line.strip()]


def load_case(db, urs, seq_len, model_length, assigned, inferred, model_id=1):
    with db.cursor() as cur:
        cur.execute("INSERT INTO rna (upi, len) VALUES (%s, %s)", (urs, seq_len))
        cur.execute(
            """INSERT INTO r2dt_models
               (id, model_name, model_source, model_length, model_basepair_count, so_term_id)
               VALUES (%s, %s, 'crw', %s, 30, 'SO:0000252')
               ON CONFLICT (id) DO NOTHING""",
            (model_id, f"model-{model_id}", model_length),
        )
        cur.execute(
            """INSERT INTO r2dt_results
               (urs, model_id, inferred_should_show, assigned_should_show,
                sequence_coverage, model_coverage, basepair_count)
               VALUES (%s, %s, %s, %s, 0.9, NULL, 25)""",
            (urs, model_id, inferred, assigned),
        )
        cur.execute(
            "INSERT INTO search_export_urs (id, urs, urs_taxid) VALUES (%s, %s, %s)",
            (model_id, urs, f"{urs}_9606"),
        )
        cur.execute(
            """INSERT INTO precompute_urs_taxid (id, precompute_urs_id, urs, urs_taxid)
               VALUES (%s, %s, %s, %s)""",
            (model_id, model_id, urs, f"{urs}_9606"),
        )


# name, sequence length, model length, assigned, inferred, should be displayed
CASES = [
    ("exact fit", 100, 100, None, True, True),
    ("just under the ratio", 124, 100, None, True, True),
    ("exactly at the ratio", 125, 100, None, True, False),
    ("force fit onto a short model", 2000, 200, None, True, False),
    ("sequence much shorter than model", 50, 3000, None, True, True),
    ("manual override shows a force fit", 2000, 200, True, False, True),
    ("manual hide beats a good ratio", 100, 100, False, True, False),
    ("null model length defers to should_show", 2000, None, None, True, True),
    ("zero model length defers to should_show", 2000, 0, None, True, True),
    ("not inferred and not assigned", 100, 100, None, False, False),
]


@pytest.mark.parametrize(
    "name,seq_len,model_length,assigned,inferred,displayed",
    CASES,
    ids=[c[0] for c in CASES],
)
@pytest.mark.parametrize(
    "query", [SEARCH_EXPORT, PRECOMPUTE], ids=["search", "precompute"]
)
def test_display_filter(
    db, query, name, seq_len, model_length, assigned, inferred, displayed
):
    load_case(db, "URS0000000001", seq_len, model_length, assigned, inferred)
    found = copy_query(db, query)
    assert bool(found) is displayed, name


def test_force_fits_are_the_only_thing_the_ratio_removes(db):
    """The ratio must not touch the sequence-shorter-than-model direction."""
    for i, (seq_len, model_length) in enumerate(
        [(50, 3000), (100, 100), (124, 100), (2000, 200)], start=1
    ):
        load_case(db, f"URS000000000{i}", seq_len, model_length, None, True, model_id=i)

    shown = {r["urs_taxid"] for r in copy_query(db, SEARCH_EXPORT)}
    assert shown == {
        "URS0000000001_9606",
        "URS0000000002_9606",
        "URS0000000003_9606",
    }


def test_both_queries_agree(db):
    for i, (seq_len, model_length, assigned, inferred, _) in enumerate(
        [(c[1], c[2], c[3], c[4], c[5]) for c in CASES], start=1
    ):
        load_case(
            db,
            f"URS00000000{i:02d}",
            seq_len,
            model_length,
            assigned,
            inferred,
            model_id=i,
        )

    search = {r["urs_taxid"] for r in copy_query(db, SEARCH_EXPORT)}
    precompute = {r["urs_taxid"] for r in copy_query(db, PRECOMPUTE)}
    assert search == precompute


def test_should_show_update_writes_back_to_r2dt_results(db):
    """The AFTER LOAD DO block must update r2dt_results, not the load table.

    It used to have the UPDATE the other way round and named a column that does
    not exist on r2dt_results: valid SQL, silent no-op, so the classifier's
    output never reached the database.
    """
    load_case(db, "URS0000000001", 100, 100, None, True, model_id=1)
    load_case(db, "URS0000000002", 100, 100, None, True, model_id=2)

    with db.cursor() as cur:
        for statement in ctl_blocks("BEFORE LOAD DO"):
            cur.execute(statement)
        cur.execute("""INSERT INTO load_secondary_should_show (urs, should_show)
               VALUES ('URS0000000001', false), ('URS0000000002', true)""")

        update = ctl_blocks("AFTER LOAD DO")[0]
        cur.execute(update)

        cur.execute("SELECT urs, inferred_should_show FROM r2dt_results ORDER BY urs")
        assert cur.fetchall() == [
            ("URS0000000001", False),
            ("URS0000000002", True),
        ]


def test_should_show_update_leaves_manual_assignments_alone(db):
    load_case(db, "URS0000000001", 100, 100, True, True, model_id=1)

    with db.cursor() as cur:
        for statement in ctl_blocks("BEFORE LOAD DO"):
            cur.execute(statement)
        cur.execute(
            "INSERT INTO load_secondary_should_show (urs, should_show) VALUES ('URS0000000001', false)"
        )
        cur.execute(ctl_blocks("AFTER LOAD DO")[0])

        cur.execute("SELECT assigned_should_show FROM r2dt_results")
        assert cur.fetchone() == (True,)

    # the manual assignment still wins in the display queries
    assert len(copy_query(db, SEARCH_EXPORT)) == 1
