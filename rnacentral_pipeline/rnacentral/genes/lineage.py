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

"""
Gene lineage planning: decide how the genes of one release relate to the genes of
the previous one, and produce the gene table, history events and member deltas.

This module is pure (polars in, polars out, no database). `plan` takes the stored
state of a taxon plus the new clusters of regions and returns:

* genes:   the next state, every gene of the taxon, active or not
* events:  one row per change (created, changed, reactivated, deactivated,
           split_from, merged_from). Genes that continue unchanged write nothing
* member_deltas: added/removed region_names per gene, keyed on
           (assembly_id, region_name) because locus ids are not stable

Identity is decided by member overlap. Every (old gene, new cluster) pair that
shares members is scored, and one greedy pass over the globally sorted pairs
decides who inherits whose identity, so the result does not depend on input order.
Extent overlap is only a fallback for old genes whose members have all gone.

`genes` has the same schema as `state`, so the output of one release can be fed
straight in as the state of the next (this is how historical GFFs are replayed).
"""

import hashlib
import uuid
from typing import NamedTuple

import polars as pl

STATE_SCHEMA = {
    "internal_name": pl.Utf8,
    "stable_id": pl.Utf8,
    "version": pl.Int64,
    "first_release": pl.Int64,
    "last_release": pl.Int64,
    "is_active": pl.Boolean,
    "assembly_id": pl.Utf8,
    "chromosome": pl.Utf8,
    "start": pl.Int64,
    "stop": pl.Int64,
    "strand": pl.Utf8,
    "members": pl.List(pl.Utf8),
    "last_members": pl.List(pl.Utf8),
}

CLUSTER_COLUMNS = [
    "cluster_id",
    "region_name",
    "assembly_id",
    "chromosome",
    "strand",
    "region_start",
    "region_stop",
]

EVENT_SCHEMA = {
    "internal_name": pl.Utf8,
    "release": pl.Int64,
    "event_type": pl.Utf8,
    "related_internal_name": pl.Utf8,
    "shared_members": pl.Int64,
    "public_name": pl.Utf8,
    "assembly_id": pl.Utf8,
    "chromosome": pl.Utf8,
    "start": pl.Int64,
    "stop": pl.Int64,
    "strand": pl.Utf8,
    "member_count": pl.Int64,
}

DELTA_SCHEMA = {
    "internal_name": pl.Utf8,
    "release": pl.Int64,
    "change": pl.Utf8,
    "assembly_id": pl.Utf8,
    "region_name": pl.Utf8,
    "urs_taxid": pl.Utf8,
}

NAMESPACE = uuid.UUID("5b0e3f0c-7a5e-4c43-a7a4-1d8a0a6f6c11")

OLD_COLUMNS = [f"o_{c}" for c in STATE_SCHEMA]


class Plan(NamedTuple):
    genes: pl.DataFrame
    events: pl.DataFrame
    member_deltas: pl.DataFrame
    stats: dict


def empty_state():
    return pl.DataFrame(schema=STATE_SCHEMA)


def _extent_overlap(a_start, a_stop, b_start, b_stop):
    """Reciprocal overlap as a polars expression: min(overlap/len_a, overlap/len_b)."""
    overlap = (
        pl.min_horizontal(a_stop, b_stop) - pl.max_horizontal(a_start, b_start)
    ).clip(lower_bound=0)
    len_a = (a_stop - a_start).clip(lower_bound=1)
    len_b = (b_stop - b_start).clip(lower_bound=1)
    return pl.min_horizontal(overlap / len_a, overlap / len_b)


def _greedy(pairs):
    """Walk pre-sorted (old, cluster) pairs, giving each old gene and cluster one partner."""
    old_taken, cluster_taken, heirs = set(), set(), {}
    for old, cluster in pairs:
        if old in old_taken or cluster in cluster_taken:
            continue
        old_taken.add(old)
        cluster_taken.add(cluster)
        heirs[cluster] = old
    return heirs


def build_clusters(clusters):
    """
    Collapse the region-level cluster table into one row per cluster.

    A cluster must sit on a single (assembly, chromosome, strand) and a region may be
    in only one cluster; anything else is a bug upstream, so raise rather than guess.
    """
    missing = [c for c in CLUSTER_COLUMNS if c not in clusters.columns]
    if missing:
        raise ValueError(f"clusters is missing columns: {missing}")

    dup = clusters.group_by("assembly_id", "region_name").agg(
        pl.col("cluster_id").n_unique().alias("n")
    )
    if dup.filter(pl.col("n") > 1).height:
        raise ValueError("a region appears in more than one cluster")

    out = (
        clusters.unique(subset=["assembly_id", "region_name"])
        .group_by("cluster_id")
        .agg(
            pl.col("assembly_id").first(),
            pl.col("chromosome").first(),
            pl.col("strand").first(),
            pl.col("region_start").min().alias("start"),
            pl.col("region_stop").max().alias("stop"),
            pl.col("region_name").sort().alias("members"),
            pl.col("assembly_id").n_unique().alias("_n_asm"),
            pl.col("chromosome").n_unique().alias("_n_chr"),
            pl.col("strand").n_unique().alias("_n_strand"),
        )
    )
    if out.filter(
        (pl.col("_n_asm") + pl.col("_n_chr") + pl.col("_n_strand")) > 3
    ).height:
        raise ValueError("a cluster spans more than one assembly, chromosome or strand")

    return (
        out.drop("_n_asm", "_n_chr", "_n_strand")
        .with_columns(
            (pl.col("assembly_id") + "|" + pl.col("members").list.first()).alias(
                "cluster_key"
            ),
            pl.col("start").cast(pl.Int64),
            pl.col("stop").cast(pl.Int64),
        )
        .drop("cluster_id")
        .sort("cluster_key")
    )


def _score_edges(state, cl):
    """One row per (old gene, cluster) pair sharing members, with the sort-key columns."""
    old = state.select(
        "internal_name",
        "assembly_id",
        pl.col("last_members").alias("region_name"),
    ).explode("region_name")
    new = cl.select(
        "cluster_key", "assembly_id", pl.col("members").alias("region_name")
    )
    new = new.explode("region_name")

    edges = (
        old.drop_nulls("region_name")
        .join(new, on=["assembly_id", "region_name"])
        .group_by("internal_name", "cluster_key")
        .agg(pl.len().alias("s"))
    )
    surv = edges.group_by("internal_name").agg(pl.col("s").sum().alias("old_surv"))
    return (
        edges.join(surv, on="internal_name")
        .join(
            state.select(
                "internal_name",
                "is_active",
                "first_release",
                pl.col("start").alias("o_start"),
                pl.col("stop").alias("o_stop"),
            ),
            on="internal_name",
        )
        .join(
            cl.select(
                "cluster_key",
                pl.col("start").alias("c_start"),
                pl.col("stop").alias("c_stop"),
                pl.col("members").list.len().alias("new_n"),
            ),
            on="cluster_key",
        )
        .with_columns(
            (pl.col("s") / (pl.col("old_surv") + pl.col("new_n") - pl.col("s"))).alias(
                "jaccard"
            ),
            _extent_overlap(
                pl.col("o_start"), pl.col("o_stop"), pl.col("c_start"), pl.col("c_stop")
            ).alias("extent"),
            (pl.col("s") / pl.col("old_surv")).alias("frac_old"),
        )
    )


def _find_heirs(state, cl, edges, min_inherit, min_extent):
    eligible = edges.filter(pl.col("frac_old") >= min_inherit).sort(
        [
            "s",
            "jaccard",
            "extent",
            "is_active",
            "first_release",
            "internal_name",
            "cluster_key",
        ],
        descending=[True, True, True, True, False, False, False],
    )
    heirs = _greedy(zip(eligible["internal_name"], eligible["cluster_key"]))

    ## Extent fallback: only active old genes that share no members with any cluster
    ## (e.g. remapped, so every region_name changed) and clusters that share none with
    ## any old gene. An inactive gene is never revived on coordinates alone.
    has_edge_old = set(edges["internal_name"])
    has_edge_new = set(edges["cluster_key"])
    left_old = state.filter(
        pl.col("is_active") & ~pl.col("internal_name").is_in(has_edge_old)
    ).select("internal_name", "assembly_id", "chromosome", "strand", "start", "stop")
    left_new = cl.filter(~pl.col("cluster_key").is_in(has_edge_new)).select(
        "cluster_key",
        "assembly_id",
        "chromosome",
        "strand",
        pl.col("start").alias("c_start"),
        pl.col("stop").alias("c_stop"),
    )
    if left_old.height and left_new.height:
        # ponytail: quadratic per (assembly, chromosome, strand) over leftovers only;
        # switch to a sorted interval sweep if the leftovers get large.
        fb = (
            left_old.join(left_new, on=["assembly_id", "chromosome", "strand"])
            .with_columns(
                _extent_overlap(
                    pl.col("start"), pl.col("stop"), pl.col("c_start"), pl.col("c_stop")
                ).alias("extent")
            )
            .filter(pl.col("extent") >= min_extent)
            .sort(
                ["extent", "internal_name", "cluster_key"],
                descending=[True, False, False],
            )
        )
        heirs.update(_greedy(zip(fb["internal_name"], fb["cluster_key"])))
    return heirs


def _mint_stable_id(assembly, chromosome, strand, start, stop, prefix, taken):
    salt = 0
    while True:
        text = f"{assembly}::{chromosome}::{strand}::{start}::{stop}::{salt}"
        digest = int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest(), "big")
        stable_id = f"RNAC{prefix}G{digest % 10**11:011d}"
        if stable_id not in taken:
            taken.add(stable_id)
            return stable_id
        salt += 1


def _mint_identities(created, release, prefix, taken):
    """Deterministic internal_name and stable_id for clusters without an heir."""
    if created.is_empty():
        return pl.DataFrame(
            schema={"cluster_key": pl.Utf8, "n_internal": pl.Utf8, "n_stable": pl.Utf8}
        )
    if prefix is None:
        raise ValueError(
            "no Ensembl stable prefix for this taxon but new genes are needed"
        )
    rows = []
    created = created.sort(
        ["assembly_id", "chromosome", "start", "stop", "strand", "cluster_key"]
    )
    for r in created.iter_rows(named=True):
        members_hash = hashlib.sha256(
            "\n".join(r["members"]).encode("utf-8")
        ).hexdigest()
        rows.append(
            {
                "cluster_key": r["cluster_key"],
                "n_internal": str(
                    uuid.uuid5(
                        NAMESPACE, f"{r['assembly_id']}|{release}|{members_hash}"
                    )
                ),
                "n_stable": _mint_stable_id(
                    r["assembly_id"],
                    r["chromosome"],
                    r["strand"],
                    r["start"],
                    r["stop"],
                    prefix,
                    taken,
                ),
            }
        )
    return pl.DataFrame(
        rows,
        schema={"cluster_key": pl.Utf8, "n_internal": pl.Utf8, "n_stable": pl.Utf8},
    )


def _event_frame(genes, release, event_type, count_col="members"):
    return genes.select(
        pl.col("internal_name"),
        pl.lit(release, dtype=pl.Int64).alias("release"),
        pl.lit(event_type).alias("event_type"),
        pl.lit(None, dtype=pl.Utf8).alias("related_internal_name"),
        pl.lit(None, dtype=pl.Int64).alias("shared_members"),
        (pl.col("stable_id") + "." + pl.col("version").cast(pl.Utf8)).alias(
            "public_name"
        ),
        "assembly_id",
        "chromosome",
        "start",
        "stop",
        "strand",
        pl.col(count_col).list.len().cast(pl.Int64).alias("member_count"),
    )


def _relation_events(edges, heirs, active_old, new_genes, release, min_relation):
    """split_from / merged_from rows for non-heir edges that carry enough members."""
    heir_pairs = pl.DataFrame(
        {"cluster_key": list(heirs.keys()), "internal_name": list(heirs.values())},
        schema={"cluster_key": pl.Utf8, "internal_name": pl.Utf8},
    )
    rel = edges.filter(pl.col("is_active") & (pl.col("frac_old") >= min_relation))
    if rel.is_empty():
        return pl.DataFrame(schema=EVENT_SCHEMA)
    rel = rel.with_columns(
        pl.len().over("internal_name").alias("n_targets"),
        pl.len().over("cluster_key").alias("n_sources"),
    )
    non_heir = rel.join(heir_pairs, on=["cluster_key", "internal_name"], how="anti")
    splits = non_heir.filter(pl.col("n_targets") >= 2).with_columns(
        pl.lit("split_from").alias("event_type")
    )
    merges = non_heir.filter(pl.col("n_sources") >= 2).with_columns(
        pl.lit("merged_from").alias("event_type")
    )
    both = pl.concat([splits, merges])
    if both.is_empty():
        return pl.DataFrame(schema=EVENT_SCHEMA)

    # the event belongs to the new gene (the cluster's gene), related to the old one
    return (
        both.join(
            new_genes.select(
                "cluster_key",
                "internal_name",
                "stable_id",
                "version",
                "assembly_id",
                "chromosome",
                "start",
                "stop",
                "strand",
                "members",
            ),
            on="cluster_key",
            suffix="_new",
        )
        .select(
            pl.col("internal_name_new").alias("internal_name"),
            pl.lit(release, dtype=pl.Int64).alias("release"),
            "event_type",
            pl.col("internal_name").alias("related_internal_name"),
            pl.col("s").cast(pl.Int64).alias("shared_members"),
            (pl.col("stable_id") + "." + pl.col("version").cast(pl.Utf8)).alias(
                "public_name"
            ),
            "assembly_id",
            "chromosome",
            "start",
            "stop",
            "strand",
            pl.col("members").list.len().cast(pl.Int64).alias("member_count"),
        )
        .unique(subset=["internal_name", "event_type", "related_internal_name"])
    )


def plan(
    state,
    clusters,
    release,
    prefix=None,
    *,
    min_inherit=0.5,
    min_relation=0.2,
    min_extent=0.9,
    allow_empty=False,
    max_deactivate_fraction=0.5,
):
    """
    Plan the transition of one taxon from `state` to the genes described by `clusters`.

    state    - STATE_SCHEMA frame: every stored gene of the taxon, active or not
    clusters - region-level frame with CLUSTER_COLUMNS
    release  - the release being produced; must not be older than anything in state
    prefix   - Ensembl stable prefix of the taxon, only needed to mint new genes
    """
    state = state.cast(STATE_SCHEMA) if state.height else empty_state()
    if state["internal_name"].n_unique() != state.height:
        raise ValueError("duplicate internal_name in state")
    if state.height and release < state["last_release"].max():
        raise ValueError(f"release {release} is older than the stored state")

    active_old = state.filter(pl.col("is_active"))
    if clusters.is_empty():
        if active_old.height and not allow_empty:
            raise ValueError("no clusters supplied but the taxon has active genes")
        cl = pl.DataFrame(
            schema={
                "assembly_id": pl.Utf8,
                "chromosome": pl.Utf8,
                "strand": pl.Utf8,
                "start": pl.Int64,
                "stop": pl.Int64,
                "members": pl.List(pl.Utf8),
                "cluster_key": pl.Utf8,
            }
        )
    else:
        cl = build_clusters(clusters)

    edges = _score_edges(state, cl) if state.height and cl.height else None
    if edges is None:
        edges = _score_edges(empty_state(), cl.clear()).clear()
    heirs = _find_heirs(state, cl, edges, min_inherit, min_extent) if cl.height else {}

    ## One row per cluster, with the heir's old state alongside
    heir_frame = pl.DataFrame(
        {"cluster_key": list(heirs.keys()), "heir": list(heirs.values())},
        schema={"cluster_key": pl.Utf8, "heir": pl.Utf8},
    )
    old = state.rename({c: f"o_{c}" for c in STATE_SCHEMA})
    rows = cl.join(heir_frame, on="cluster_key", how="left").join(
        old, left_on="heir", right_on="o_internal_name", how="left"
    )

    created = rows.filter(pl.col("heir").is_null())
    taken = set(state["stable_id"].to_list())
    minted = _mint_identities(created, release, prefix, taken)
    rows = rows.join(minted, on="cluster_key", how="left")

    same_active = pl.col("o_members").list.sort() == pl.col("members")
    same_last = pl.col("o_last_members").list.sort() == pl.col("members")
    status = (
        pl.when(pl.col("heir").is_null())
        .then(pl.lit("created"))
        .when(pl.col("o_is_active") & same_active)
        .then(pl.lit("continued"))
        .when(pl.col("o_is_active"))
        .then(pl.lit("changed"))
        .otherwise(pl.lit("reactivated"))
    )
    version = (
        pl.when(pl.col("status") == "created")
        .then(1)
        .when(pl.col("status") == "continued")
        .then(pl.col("o_version"))
        .when(pl.col("status") == "changed")
        .then(pl.col("o_version") + 1)
        .when(same_last)
        .then(pl.col("o_version"))
        .otherwise(pl.col("o_version") + 1)
    )
    new_genes = (
        rows.with_columns(status.alias("status"))
        .with_columns(version.cast(pl.Int64).alias("version"))
        .select(
            "cluster_key",
            pl.coalesce("heir", "n_internal").alias("internal_name"),
            pl.coalesce("o_stable_id", "n_stable").alias("stable_id"),
            "version",
            pl.coalesce("o_first_release", pl.lit(release, dtype=pl.Int64)).alias(
                "first_release"
            ),
            pl.lit(release, dtype=pl.Int64).alias("last_release"),
            pl.lit(True).alias("is_active"),
            "assembly_id",
            "chromosome",
            "start",
            "stop",
            "strand",
            "members",
            pl.col("members").alias("last_members"),
            "status",
        )
    )

    heir_names = set(heirs.values())
    deactivated = state.filter(
        pl.col("is_active") & ~pl.col("internal_name").is_in(heir_names)
    )
    if (
        active_old.height
        and deactivated.height / active_old.height > max_deactivate_fraction
        and not (allow_empty and cl.is_empty())
    ):
        raise ValueError(
            f"{deactivated.height} of {active_old.height} active genes would be "
            f"deactivated, over the {max_deactivate_fraction} limit"
        )
    gone = deactivated.with_columns(
        pl.lit(False).alias("is_active"),
        pl.col("members").alias("last_members"),
        pl.lit([], dtype=pl.List(pl.Utf8)).alias("members"),
        pl.lit("deactivated").alias("status"),
    )
    still_inactive = state.filter(
        ~pl.col("is_active") & ~pl.col("internal_name").is_in(heir_names)
    ).with_columns(pl.lit("inactive").alias("status"))

    genes = pl.concat(
        [
            new_genes.drop("cluster_key").select(list(STATE_SCHEMA) + ["status"]),
            gone.select(list(STATE_SCHEMA) + ["status"]),
            still_inactive.select(list(STATE_SCHEMA) + ["status"]),
        ]
    ).sort("internal_name")
    genes = genes.with_columns(
        (pl.col("stable_id") + "." + pl.col("version").cast(pl.Utf8)).alias(
            "public_name"
        )
    )

    ## Events
    snap = new_genes.drop("cluster_key")
    events = pl.concat(
        [
            _event_frame(
                snap.filter(pl.col("status") == "created"), release, "created"
            ),
            _event_frame(
                snap.filter(pl.col("status") == "changed"), release, "changed"
            ),
            _event_frame(
                snap.filter(pl.col("status") == "reactivated"), release, "reactivated"
            ),
            _event_frame(gone, release, "deactivated", count_col="last_members"),
            _relation_events(
                edges, heirs, active_old, new_genes, release, min_relation
            ),
        ]
    ).sort(["internal_name", "event_type", "related_internal_name"], nulls_last=True)

    ## Member deltas: compare every active old gene's members with the new gene's
    p = (
        active_old.select(
            "internal_name", "assembly_id", pl.col("members").alias("region_name")
        )
        .explode("region_name")
        .drop_nulls("region_name")
    )
    n = (
        snap.select(
            "internal_name", "assembly_id", pl.col("members").alias("region_name")
        )
        .explode("region_name")
        .drop_nulls("region_name")
    )
    key = ["internal_name", "region_name"]
    deltas = (
        pl.concat(
            [
                n.join(p, on=key, how="anti").with_columns(
                    pl.lit("added").alias("change")
                ),
                p.join(n, on=key, how="anti").with_columns(
                    pl.lit("removed").alias("change")
                ),
            ]
        )
        .select(
            "internal_name",
            pl.lit(release, dtype=pl.Int64).alias("release"),
            "change",
            "assembly_id",
            "region_name",
            pl.col("region_name").str.split("@").list.first().alias("urs_taxid"),
        )
        .sort(["internal_name", "change", "region_name"])
    )

    stats = {
        status: int((genes["status"] == status).sum())
        for status in ("created", "continued", "changed", "reactivated", "deactivated")
    }
    return Plan(genes, events, deltas, stats)


def clusters_from_assignments(transcripts, assignments, assembly_id=None):
    """
    Turn gene assignments (region_name, gene_id) from the GNN into the cluster frame
    that plan() takes, using the transcripts frame produced by `genes utils convert`.

    This is the one place strand changes from 1/-1 to "+"/"-". Transcripts the GNN did
    not assign (single-member clusters are never genes) are simply not in the result.
    Pass assembly_id when transcripts holds several assemblies; the GNN keys its nodes
    on region_name alone, so it must have been run for one assembly at a time.
    """
    if assembly_id is not None:
        transcripts = transcripts.filter(pl.col("assembly_id") == assembly_id)
    transcripts = transcripts.select(
        "assembly_id",
        "chromosome",
        "region_name",
        "region_start",
        "region_stop",
        "strand",
    ).unique(subset=["assembly_id", "region_name"])
    if transcripts["region_name"].n_unique() != transcripts.height:
        raise ValueError(
            "a region_name is on more than one assembly; pass assembly_id and run the "
            "GNN once per assembly"
        )
    assignments = assignments.select("region_name", "gene_id").unique()
    joined = assignments.join(transcripts, on="region_name", how="left")
    unknown = joined.filter(pl.col("assembly_id").is_null())
    if unknown.height:
        raise ValueError(
            f"{unknown.height} assigned regions are not in the transcripts, "
            f"e.g. {unknown['region_name'][0]}"
        )
    return joined.select(
        (pl.col("assembly_id") + "|" + pl.col("gene_id")).alias("cluster_id"),
        "region_name",
        "assembly_id",
        "chromosome",
        pl.when(pl.col("strand") == 1)
        .then(pl.lit("+"))
        .otherwise(pl.lit("-"))
        .alias("strand"),
        pl.col("region_start").cast(pl.Int64),
        pl.col("region_stop").cast(pl.Int64),
    )


def as_state(genes):
    """The genes frame of one plan, ready to be the state of the next release."""
    return genes.select(list(STATE_SCHEMA))
