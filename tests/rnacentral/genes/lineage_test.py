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

import random

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.genes import lineage as L

PREFIX = "TEST"


def region(i, start=None, stop=None, strand="+"):
    start = 1000 * i if start is None else start
    stop = start + 100 if stop is None else stop
    return f"URS{i:010d}_9606@1/{start}-{stop}:{strand}"


def clusters(*groups, assembly="A1", strand="+"):
    """groups: lists of region indices (or region_name strings); one cluster each."""
    rows = []
    for gi, group in enumerate(groups):
        for m in group:
            name = m if isinstance(m, str) else region(m, strand=strand)
            start, stop = [int(x) for x in name.split("/")[1].split(":")[0].split("-")]
            rows.append(
                dict(
                    cluster_id=f"c{gi}",
                    region_name=name,
                    assembly_id=assembly,
                    chromosome="1",
                    strand=strand,
                    region_start=start,
                    region_stop=stop,
                )
            )
    return (
        pl.DataFrame(
            rows, schema_overrides={"region_start": pl.Int64, "region_stop": pl.Int64}
        )
        if rows
        else pl.DataFrame(schema={c: pl.Utf8 for c in L.CLUSTER_COLUMNS})
    )


def run(state, groups, release, **kw):
    return L.plan(state, clusters(*groups, **kw.pop("cl", {})), release, PREFIX, **kw)


def chain(releases, start=1):
    """releases: list of group lists. Returns the plans, feeding each into the next."""
    state, plans = L.empty_state(), []
    for offset, groups in enumerate(releases):
        p = run(state, groups, start + offset)
        plans.append(p)
        state = L.as_state(p.genes)
    return plans


def active(genes):
    return genes.filter(pl.col("is_active"))


def members_of(genes, internal_name):
    return set(genes.filter(pl.col("internal_name") == internal_name)["members"][0])


def regions(*ids):
    return {region(i) for i in ids}


def test_first_release_creates_everything():
    p = run(L.empty_state(), [[1, 2, 3], [10, 11]], 1)
    assert p.stats["created"] == 2
    assert set(p.events["event_type"]) == {"created"}
    assert p.member_deltas["change"].unique().to_list() == ["added"]
    assert p.member_deltas.height == 5
    g = active(p.genes)
    assert (g["version"] == 1).all() and (g["first_release"] == 1).all()
    assert g["public_name"].str.contains(r"^RNACTESTG\d{11}\.1$").all()


def test_rerun_same_release_is_a_noop():
    p1 = run(L.empty_state(), [[1, 2, 3], [10, 11]], 1)
    p2 = run(L.as_state(p1.genes), [[1, 2, 3], [10, 11]], 1)
    assert p2.events.is_empty() and p2.member_deltas.is_empty()
    assert p2.genes.equals(p1.genes.with_columns(pl.lit("continued").alias("status")))


def test_unchanged_next_release_writes_nothing_but_last_release():
    p1, p2 = chain([[[1, 2, 3]], [[1, 2, 3]]])
    assert p2.events.is_empty() and p2.member_deltas.is_empty()
    g = p2.genes
    assert g["last_release"].to_list() == [2] and g["first_release"].to_list() == [1]
    assert g["version"].to_list() == [1]


def test_member_change_keeps_identity_and_bumps_version():
    p1, p2 = chain([[[1, 2, 3, 4]], [[1, 2, 3, 5]]])
    assert p2.genes["internal_name"].to_list() == p1.genes["internal_name"].to_list()
    assert p2.genes["version"].to_list() == [2]
    assert p2.events["event_type"].to_list() == ["changed"]
    assert sorted(p2.member_deltas.select("change", "region_name").rows()) == [
        ("added", region(5)),
        ("removed", region(4)),
    ]


def test_extent_is_the_new_clusters_own_and_does_not_ratchet():
    p1, p2, p3 = chain([[[1, 2]], [[1, 2, 3]], [[1, 2]]])
    assert p2.genes["stop"][0] > p1.genes["stop"][0]
    assert p3.genes["stop"][0] == p1.genes["stop"][0]


def test_stable_id_survives_member_changes():
    p1, p2 = chain([[[1, 2, 3, 4]], [[1, 2, 3, 5]]])
    assert p1.genes["stable_id"].to_list() == p2.genes["stable_id"].to_list()


def test_split_50_50_is_deterministic_and_records_the_split():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3, 4]], 1).genes)
    a = run(state, [[1, 2], [3, 4]], 2)
    b = run(state, [[3, 4], [1, 2]], 2)
    assert (
        a.genes.drop("internal_name", "stable_id", "public_name")
        .sort("members")
        .equals(
            b.genes.drop("internal_name", "stable_id", "public_name").sort("members")
        )
    )
    assert a.stats["created"] == 1 and a.stats["changed"] == 1
    split = a.events.filter(pl.col("event_type") == "split_from")
    assert split.height == 1
    assert split["related_internal_name"][0] == state["internal_name"][0]
    assert split["shared_members"][0] == 2


def test_split_heir_is_the_larger_side():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3, 4, 5]], 1).genes)
    p = run(state, [[1, 2, 3], [4, 5]], 2)
    heir = p.genes.filter(pl.col("internal_name") == state["internal_name"][0])
    assert members_of(heir, state["internal_name"][0]) == regions(1, 2, 3)


def test_merge_two_into_one():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3], [4, 5, 6]], 1).genes)
    p = run(state, [[1, 2, 3, 4, 5, 6]], 2)
    assert p.stats["deactivated"] == 1 and p.stats["changed"] == 1
    merged = p.events.filter(pl.col("event_type") == "merged_from")
    assert merged.height == 1
    assert merged["related_internal_name"][0] in set(state["internal_name"])
    assert "deactivated" in set(p.events["event_type"])


def test_many_to_many_gives_split_and_merge_rows():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3, 4], [5, 6, 7, 8]], 1).genes)
    p = run(state, [[1, 2, 3, 5, 6], [4, 7, 8]], 2)
    types = set(p.events["event_type"])
    assert {"split_from", "merged_from"} <= types


def test_singleton_move_is_recorded_in_member_deltas():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3], [4, 5, 6]], 1).genes)
    p = run(state, [[1, 2], [3, 4, 5, 6]], 2)
    moved = p.member_deltas.filter(pl.col("region_name") == region(3))
    assert sorted(moved["change"].to_list()) == ["added", "removed"]
    assert moved["internal_name"].n_unique() == 2


def test_below_inherit_threshold_is_a_new_gene():
    state = L.as_state(run(L.empty_state(), [[1, 2, 3, 4, 5, 6]], 1).genes)
    p = run(state, [[1, 7, 8, 9, 10, 11], [2, 3, 4, 5, 6]], 2)
    ## the 5-of-6 cluster inherits; the 1-of-6 cluster does not
    heir = p.genes.filter(pl.col("internal_name") == state["internal_name"][0])
    assert members_of(heir, state["internal_name"][0]) == regions(2, 3, 4, 5, 6)
    assert p.stats["created"] == 1


def test_gene_that_grows_still_inherits():
    state = L.as_state(run(L.empty_state(), [[1, 2]], 1).genes)
    p = run(state, [list(range(1, 11))], 2)
    assert p.stats["changed"] == 1 and p.stats["created"] == 0


def test_extent_fallback_when_every_region_was_remapped():
    old = [region(1, 1000, 1100), region(2, 1200, 1300)]
    new = [region(1, 1001, 1101), region(2, 1201, 1301)]
    state = L.as_state(run(L.empty_state(), [old], 1).genes)
    p = run(state, [new], 2)
    assert p.stats["created"] == 0 and p.stats["changed"] == 1
    assert p.genes["internal_name"].to_list() == state["internal_name"].to_list()


def test_extent_fallback_never_revives_an_inactive_gene():
    s1 = L.as_state(run(L.empty_state(), [[1, 2]], 1).genes)
    s2 = L.as_state(run(s1, [], 2, allow_empty=True).genes)
    assert not s2["is_active"][0]
    remapped = [region(1, 1001, 1101), region(2, 1201, 1301)]
    p = run(s2, [remapped], 3)
    assert p.stats["created"] == 1 and p.stats["reactivated"] == 0


def test_deactivated_gene_reactivates_with_same_identity():
    s1 = L.as_state(run(L.empty_state(), [[1, 2, 3]], 1).genes)
    s2 = L.as_state(run(s1, [[10, 11]], 2, max_deactivate_fraction=1.0).genes)
    assert (
        s2.filter(pl.col("internal_name") == s1["internal_name"][0])["is_active"][0]
        is False
    )
    p = run(s2, [[10, 11], [1, 2, 3]], 3)
    assert p.stats["reactivated"] == 1
    assert s1["internal_name"][0] in set(active(p.genes)["internal_name"])
    ## set unchanged since deactivation, so the version is not bumped
    row = p.genes.filter(pl.col("internal_name") == s1["internal_name"][0])
    assert row["version"][0] == s1["version"][0]


def test_deactivation_logs_removed_for_every_member():
    s1 = L.as_state(run(L.empty_state(), [[1, 2, 3], [10]], 1).genes)
    p = run(s1, [[10]], 2, max_deactivate_fraction=0.9)
    removed = p.member_deltas.filter(pl.col("change") == "removed")
    assert set(removed["region_name"]) == regions(1, 2, 3)


def test_active_beats_inactive_on_a_tie():
    s1 = L.as_state(run(L.empty_state(), [[1, 2]], 1).genes)
    ## make a second, inactive gene that also holds the same regions
    dead = s1.with_columns(
        pl.lit("dead").alias("internal_name"),
        pl.lit("RNACTESTG00000000001").alias("stable_id"),
        pl.lit(False).alias("is_active"),
        pl.lit([], dtype=pl.List(pl.Utf8)).alias("members"),
    )
    p = run(pl.concat([s1, dead]), [[1, 2]], 2)
    assert s1["internal_name"][0] in set(active(p.genes)["internal_name"])
    assert "dead" not in set(active(p.genes)["internal_name"])


def test_memberless_old_gene_is_kept_not_duplicated():
    s1 = L.as_state(run(L.empty_state(), [[1, 2]], 1).genes)
    memberless = s1.with_columns(
        pl.lit([], dtype=pl.List(pl.Utf8)).alias("members"),
        pl.lit([], dtype=pl.List(pl.Utf8)).alias("last_members"),
    )
    ## members all gone from the database but the extent still matches
    p = run(memberless, [[1, 2]], 2)
    assert p.genes.height == 1 and p.stats["changed"] == 1


def test_empty_clusters_need_allow_empty():
    s1 = L.as_state(run(L.empty_state(), [[1, 2]], 1).genes)
    with pytest.raises(ValueError, match="no clusters"):
        run(s1, [], 2)
    p = run(s1, [], 2, allow_empty=True)
    assert p.stats["deactivated"] == 1


def test_deactivate_fraction_guard():
    s1 = L.as_state(run(L.empty_state(), [[1], [2], [3], [4]], 1).genes)
    with pytest.raises(ValueError, match="limit"):
        run(s1, [[1]], 2)


def test_assembly_switch_does_not_trip_the_guard_but_mass_loss_on_a_kept_assembly_does():
    old = L.as_state(
        run(L.empty_state(), [[1, 2], [3, 4], [5, 6]], 1, cl={"assembly": "A1"}).genes
    )
    ## the whole taxon moves to a new assembly: every A1 gene retires, which is expected
    p = L.plan(old, clusters([1, 2], [3, 4], assembly="A2"), 2, PREFIX)
    assert p.stats["deactivated"] == 3 and p.stats["created"] == 2
    ## but losing most of an assembly that is still present is still refused
    with pytest.raises(ValueError, match="limit"):
        L.plan(old, clusters([1, 2], assembly="A1"), 2, PREFIX)


def test_cluster_spanning_two_strands_raises():
    cl = pl.concat(
        [
            clusters([1]),
            clusters([2], strand="-").with_columns(pl.lit("c0").alias("cluster_id")),
        ]
    )
    with pytest.raises(ValueError, match="more than one"):
        L.plan(L.empty_state(), cl, 1, PREFIX)


def test_region_in_two_clusters_raises():
    with pytest.raises(ValueError, match="more than one cluster"):
        L.plan(L.empty_state(), clusters([1, 2], [2, 3]), 1, PREFIX)


def test_same_region_name_in_two_assemblies_stays_two_members():
    cl = pl.concat(
        [
            clusters([1], assembly="A1"),
            clusters([1], assembly="A2").with_columns(pl.lit("c1").alias("cluster_id")),
        ]
    )
    p = L.plan(L.empty_state(), cl, 1, PREFIX)
    assert p.genes.height == 2
    assert set(p.genes["assembly_id"]) == {"A1", "A2"}


def test_going_backwards_in_release_raises():
    s1 = L.as_state(run(L.empty_state(), [[1, 2]], 5).genes)
    with pytest.raises(ValueError, match="older"):
        run(s1, [[1, 2]], 4)


def test_new_genes_need_a_prefix():
    with pytest.raises(ValueError, match="prefix"):
        L.plan(L.empty_state(), clusters([1, 2]), 1, None)


def test_strand_is_part_of_the_name():
    cl = pl.concat(
        [
            clusters([region(1, 100, 200, "+")], strand="+"),
            clusters([region(2, 100, 200, "-")], strand="-").with_columns(
                pl.lit("c1").alias("cluster_id")
            ),
        ]
    )
    p = L.plan(L.empty_state(), cl, 1, PREFIX)
    assert p.genes["stable_id"].n_unique() == 2


def test_name_does_not_depend_on_other_chromosomes():
    one = clusters([1, 2])
    other = clusters([5, 6]).with_columns(
        pl.lit("0").alias("chromosome"), pl.lit("c9").alias("cluster_id")
    )
    a = L.plan(L.empty_state(), one, 1, PREFIX)
    b = L.plan(L.empty_state(), pl.concat([one, other]), 1, PREFIX)
    assert set(a.genes["stable_id"]) <= set(b.genes["stable_id"])


def test_extent_uses_every_exon_not_the_first():
    name = "URS0000000001_9606@1/100-200,300-400:+"
    cl = pl.DataFrame(
        [
            dict(
                cluster_id="c",
                region_name=name,
                assembly_id="A1",
                chromosome="1",
                strand="+",
                region_start=100,
                region_stop=400,
            )
        ]
    )
    p = L.plan(L.empty_state(), cl, 1, PREFIX)
    assert (p.genes["start"][0], p.genes["stop"][0]) == (100, 400)


def test_minting_is_deterministic():
    a = run(L.empty_state(), [[1, 2], [3, 4]], 1)
    b = run(L.empty_state(), [[3, 4], [1, 2]], 1)
    assert a.genes.equals(b.genes)


## Property tests over random multi-release histories


def random_history(seed, n_releases=4, n_regions=60):
    rng = random.Random(seed)
    universe = list(range(1, n_regions + 1))
    releases = []
    for _ in range(n_releases):
        pool = rng.sample(universe, rng.randint(n_regions // 2, n_regions))
        rng.shuffle(pool)
        groups, i = [], 0
        while i < len(pool):
            size = rng.randint(1, 6)
            groups.append(pool[i : i + size])
            i += size
        releases.append(groups)
    return releases


def run_relaxed(state, groups, release):
    return L.plan(
        state, clusters(*groups), release, PREFIX, max_deactivate_fraction=1.0
    )


@pytest.mark.parametrize("seed", range(40))
def test_history_properties(seed):
    state = L.empty_state()
    replay = {}  # (internal_name, region) -> present, rebuilt from deltas
    for release, groups in enumerate(random_history(seed), start=1):
        p = run_relaxed(state, groups, release)

        ## idempotent: planning the same release again changes nothing
        again = run_relaxed(L.as_state(p.genes), groups, release)
        assert again.events.is_empty() and again.member_deltas.is_empty()

        ## order independent
        shuffled = [list(reversed(g)) for g in reversed(groups)]
        other = run_relaxed(state, shuffled, release)
        assert other.genes.equals(p.genes)
        assert other.events.equals(p.events)
        assert other.member_deltas.equals(p.member_deltas)

        ## no region sits in two active genes
        mem = (
            active(p.genes)
            .select("assembly_id", pl.col("members").alias("r"))
            .explode("r")
        )
        assert mem.height == mem.unique().height

        ## replaying the deltas gives exactly the current membership
        for row in p.member_deltas.iter_rows(named=True):
            key = (row["internal_name"], row["region_name"])
            if row["change"] == "added":
                assert key not in replay
                replay[key] = True
            else:
                assert key in replay
                del replay[key]
        expected = {
            (r["internal_name"], m)
            for r in active(p.genes).iter_rows(named=True)
            for m in r["members"]
        }
        assert set(replay) == expected

        ## every old active gene is an heir or was deactivated
        old_active = set(state.filter(pl.col("is_active"))["internal_name"])
        new_active = set(active(p.genes)["internal_name"])
        deactivated = set(
            p.genes.filter(pl.col("status") == "deactivated")["internal_name"]
        )
        assert old_active <= new_active | deactivated
        assert not (new_active & deactivated)

        state = L.as_state(p.genes)
