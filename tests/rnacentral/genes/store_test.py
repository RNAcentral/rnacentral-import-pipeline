# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.genes import lineage as L
from rnacentral_pipeline.rnacentral.genes import store
from tests.rnacentral.genes.lineage_test import PREFIX, clusters, regions, run


def gene_rows(*rows):
    return pl.DataFrame(
        rows,
        schema={
            "internal_name": pl.Utf8,
            "stable_id": pl.Utf8,
            "version": pl.Int64,
            "first_release": pl.Int64,
            "last_release": pl.Int64,
            "assembly_id": pl.Utf8,
            "chromosome": pl.Utf8,
            "start": pl.Int64,
            "stop": pl.Int64,
            "strand": pl.Utf8,
            "status_event": pl.Utf8,
            "status_release": pl.Int64,
        },
        orient="row",
    )


def member_rows(*rows):
    return pl.DataFrame(
        rows,
        schema={
            "internal_name": pl.Utf8,
            "region_name": pl.Utf8,
            "change": pl.Utf8,
            "release": pl.Int64,
        },
        orient="row",
    )


def g(name, event, release=1):
    return (name, f"RNAC{name}", 1, 1, 3, "A1", "1", 100, 200, "+", event, release)


def test_active_gene_gets_the_members_whose_latest_change_is_added():
    state = store.assemble_state(
        gene_rows(g("a", "created")),
        member_rows(
            ("a", "r1", "added", 1), ("a", "r2", "added", 2), ("a", "r3", "removed", 2)
        ),
    )
    row = state.row(0, named=True)
    assert row["is_active"] and row["members"] == ["r1", "r2"]
    assert row["last_members"] == ["r1", "r2"]


def test_inactive_gene_keeps_only_what_it_lost_when_it_was_deactivated():
    state = store.assemble_state(
        gene_rows(g("a", "deactivated", release=3)),
        member_rows(
            ("a", "r1", "removed", 3),  # lost with the gene
            ("a", "r2", "removed", 2),  # moved to another gene earlier
        ),
    )
    row = state.row(0, named=True)
    assert not row["is_active"] and row["members"] == []
    assert row["last_members"] == ["r1"]


def test_reactivated_gene_is_active_again():
    state = store.assemble_state(
        gene_rows(g("a", "reactivated", release=4)),
        member_rows(("a", "r1", "added", 4)),
    )
    assert state["is_active"].to_list() == [True]
    assert state["members"].to_list() == [["r1"]]


def test_memberless_gene_survives():
    state = store.assemble_state(gene_rows(g("a", "created")), member_rows())
    assert state.height == 1 and state["members"].to_list() == [[]]


def test_gene_without_history_is_an_error_not_a_guess():
    with pytest.raises(ValueError, match="no history"):
        store.assemble_state(gene_rows(g("a", None, None)), member_rows())


def test_gene_without_stable_id_is_an_error():
    rows = gene_rows(g("a", "created")).with_columns(
        pl.lit(None, dtype=pl.Utf8).alias("stable_id")
    )
    with pytest.raises(ValueError, match="stable_id"):
        store.assemble_state(rows, member_rows())


def test_no_genes_is_an_empty_state():
    assert store.assemble_state(gene_rows(), member_rows()).equals(L.empty_state())


def test_state_from_history_round_trips_through_the_planner():
    """History written by a plan, read back as state, plans to nothing new."""
    p = run(L.empty_state(), [[1, 2, 3], [10, 11]], 1)
    gene_rows_ = p.genes.select(
        "internal_name",
        "stable_id",
        "version",
        "first_release",
        "last_release",
        "assembly_id",
        "chromosome",
        "start",
        "stop",
        "strand",
        pl.lit("created").alias("status_event"),
        pl.lit(1, dtype=pl.Int64).alias("status_release"),
    )
    members = p.member_deltas.select(
        "internal_name", "region_name", "change", "release"
    )
    state = store.assemble_state(gene_rows_, members)
    again = L.plan(state, clusters([1, 2, 3], [10, 11]), 1, PREFIX)
    assert again.events.is_empty() and again.member_deltas.is_empty()


def test_genes_for_upsert_skips_untouched_inactive_genes_and_counts_members():
    s1 = L.as_state(run(L.empty_state(), [[1, 2], [10]], 1).genes)
    s2 = L.as_state(run(s1, [[10]], 2, max_deactivate_fraction=1.0).genes)
    p3 = run(s2, [[10]], 3)
    assert p3.genes.filter(pl.col("status") == "inactive").height == 1
    rows = store.genes_for_upsert(p3.genes, 9606)
    assert rows.height == 1 and rows["member_count"].to_list() == [1]
    assert rows["taxid"].to_list() == [9606]
    p2 = run(s1, [[10]], 2, max_deactivate_fraction=1.0)
    rows2 = store.genes_for_upsert(p2.genes, 9606)
    assert sorted(rows2["member_count"].to_list()) == [0, 1]
    assert rows2.filter(~pl.col("is_active"))["member_count"].to_list() == [0]


def test_members_for_load_is_only_active_members():
    s1 = L.as_state(run(L.empty_state(), [[1, 2], [10]], 1).genes)
    p2 = run(s1, [[10]], 2, max_deactivate_fraction=1.0)
    m = store.members_for_load(p2.genes)
    assert m["region_name"].to_list() == sorted(regions(10))


def test_sql_constants_use_only_named_parameters_consistently():
    for sql in [
        store.GENE_ROWS,
        store.MEMBER_ROWS,
        store.REF_MEMBER_ROWS,
        store.PREFIX,
        store.DELETE_STALE_MEMBERS,
    ]:
        assert "%(taxid)s" in sql
    for sql in store.EVIDENCE.values():
        assert "%(taxid)s" in sql


def test_include_inactive_writes_genes_retired_earlier():
    s1 = L.as_state(run(L.empty_state(), [[1, 2], [10]], 1).genes)
    s2 = L.as_state(run(s1, [[10]], 2, max_deactivate_fraction=1.0).genes)
    p3 = run(s2, [[10]], 3)
    assert store.genes_for_upsert(p3.genes, 9606).height == 1
    both = store.genes_for_upsert(p3.genes, 9606, include_inactive=True)
    assert both.height == 2 and both["is_active"].to_list().count(False) == 1


def stored(plans):
    """
    What fetch_state reads back after `plans` were applied in order: emulates GENE_ROWS,
    MEMBER_ROWS and REF_MEMBER_ROWS over the rows apply_plan would have written.
    """
    last = plans[-1]
    events = pl.concat([p.events for p in plans])
    status = (
        events.filter(
            pl.col("event_type").is_in(["created", "reactivated", "deactivated"])
        )
        .sort("release")
        .group_by("internal_name")
        .agg(
            pl.col("event_type").last().alias("status_event"),
            pl.col("release").last().alias("status_release"),
        )
    )
    genes = store.genes_for_upsert(last.genes, 9606, include_inactive=True)
    gene_rows_ = genes.join(status, on="internal_name", how="left").drop(
        "public_name", "member_count", "taxid", "is_active"
    )
    history = pl.concat([p.member_deltas for p in plans]).select(
        "internal_name", "region_name", "change", "release"
    )

    def latest(h):
        return h.sort("release").group_by("internal_name", "region_name").last()

    refs = gene_rows_.filter(pl.col("ref_release") < pl.col("last_release"))
    ref_rows = (
        latest(
            history.join(
                refs.select("internal_name", "ref_release"), on="internal_name"
            )
            .filter(pl.col("release") <= pl.col("ref_release"))
            .drop("ref_release")
        )
        .filter(pl.col("change") == "added")
        .select("internal_name", "region_name")
    )
    return store.assemble_state(gene_rows_, latest(history), ref_rows)


def test_reference_membership_survives_the_database():
    """A gene cut down for one release takes its cluster back when read from the database,
    exactly as it does with the planner's in-memory state (the PR #315 review finding)."""
    full = list(range(1, 101))
    p1 = run(L.empty_state(), [full], 1)
    p2 = run(L.as_state(p1.genes), [[1, 2, 3]], 2)
    name = p1.genes["internal_name"][0]

    state = stored([p1, p2])
    row = state.row(0, named=True)
    assert row["ref_release"] == 1 and len(row["ref_members"]) == 100
    memory = L.as_state(p2.genes).select(state.columns).sort("internal_name")
    assert state.sort("internal_name").equals(memory)

    p3 = run(state, [full], 3, max_deactivate_fraction=1.0)
    assert set(p3.genes.filter(pl.col("is_active"))["internal_name"]) == {name}
    assert p3.stats["created"] == 0


def test_null_ref_release_means_no_older_reference():
    state = store.assemble_state(
        gene_rows(g("a", "created")), member_rows(("a", "r1", "added", 1))
    )
    row = state.row(0, named=True)
    assert row["ref_release"] == 3 and row["ref_members"] == ["r1"]
