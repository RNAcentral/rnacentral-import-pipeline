# -*- coding: utf-8 -*-

import json

import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.genes.random_forest import data, merge


def region(urs, start, stop, chrom="1", strand="+"):
    return f"{urs}_9606@{chrom}/{start}-{stop}:{strand}"


def gene(name, members, release=26, internal=None):
    return {
        "name": name,
        "internal_name": internal or "u-" + name,
        "assembly_id": "GRCh38",
        "chromosome": "1",
        "strand": "+",
        "members": members,
        "first_release": release,
        "last_release": release,
    }


A, B, C = region("URS1", 100, 200), region("URS2", 150, 260), region("URS3", 5000, 5100)
D, E = region("URS4", 90, 3000), region("URS5", 5050, 5200)


def run(previous, following, **kwargs):
    genes, events = merge.merge(previous, following, **kwargs)
    return {g["internal_name"]: g for g in genes}, events


def test_region_helpers():
    name = "URS000_9606@KN1.2/10-20,30-45:-"
    assert merge.region_span(name) == (10, 45)
    assert merge.region_urs(name) == "URS000"


def test_growth_keeps_identity_and_recomputes_the_extent():
    out, events = run([gene("g1.1", [A, B])], [gene("n.1", [A, B, D], release=28)])
    assert list(out) == ["u-g1.1"]
    assert (out["u-g1.1"]["name"], out["u-g1.1"]["start"], out["u-g1.1"]["stop"]) == (
        "g1.2",
        90,
        3000,
    )
    assert out["u-g1.1"]["first_release"] == 26 and out["u-g1.1"]["last_release"] == 28
    assert events["changed"] == 1


def test_same_members_keeps_the_version():
    out, _ = run([gene("g1.3", [A, B])], [gene("n.1", [B, A], release=28)])
    assert out["u-g1.3"]["name"] == "g1.3"


def test_a_small_fragment_does_not_take_the_identity_of_a_large_gene():
    big = [region(f"URS{i}", i * 10, i * 10 + 5) for i in range(10)]
    out, events = run(
        [gene("big.1", big)], [gene("rest.1", big[2:]), gene("frag.1", big[:2])]
    )
    assert out["u-big.1"]["members"] == sorted(big[2:])
    assert (
        out["u-frag.1"]["history"]
        == [{"event": "split_from", "internal_name": "u-big.1", "name": "big.1"}]
        and out["u-frag.1"]["version"] == 1
    )
    assert events["new_from_split"] == 1 and events["split"] == 1


def test_an_even_three_way_split_retires_the_original():
    ms = [region(f"URS{i}", i * 10, i * 10 + 5) for i in range(9)]
    out, events = run(
        [gene("big.1", ms)], [gene(f"p{i}.1", ms[3 * i : 3 * i + 3]) for i in range(3)]
    )
    assert "u-big.1" not in out and events["retired"] == 1
    assert all(
        [h["name"] for h in g["history"]] == ["big.1"]
        and g["history"][0]["event"] == "split_from"
        for g in out.values()
    )


def test_merge_history_names_the_absorbed_gene():
    prev = [gene("g1.1", [A, B, D]), gene("g2.1", [C])]
    out, events = run(prev, [gene("n.1", [A, B, D, C, E], release=28)])
    assert list(out) == ["u-g1.1"] and out["u-g1.1"]["history"] == [
        {"event": "merged_from", "internal_name": "u-g2.1", "name": "g2.1"}
    ]
    assert events["retired"] == 1  # g2 has no heir of its own


def test_retired_urs_shrink_and_retire_genes():
    prev = [gene("g1.1", [A, B]), gene("g2.1", [C])]
    out, events = run(
        prev, [gene("n1.1", [A], release=28)], retired_urs={"URS2", "URS3"}
    )
    assert (
        list(out) == ["u-g1.1"]
        and out["u-g1.1"]["name"] == "g1.2"
        and events["retired"] == 1
    )


def test_genes_on_different_assemblies_never_match():
    other = {**gene("n.1", [A, B], release=28), "assembly_id": "GRCh37"}
    out, _ = run([gene("g1.1", [A, B])], [other])
    assert list(out) == ["u-n.1"]


def test_merge_genes_end_to_end(tmp_path):
    prev, nxt, inactive = (
        tmp_path / "p.json",
        tmp_path / "n.json",
        tmp_path / "inactive.csv",
    )
    pl.DataFrame([gene("g1.1", [A, B, C]) | {"version": 1}]).write_json(prev)
    following = gene("n.1", [A, B], release=28) | {"start": 1, "stop": 2, "version": 1}
    pl.DataFrame([following]).write_json(nxt)
    inactive.write_text("URS3\n")
    out = data.merge_genes(str(prev), str(nxt), None, str(inactive), 26, 28)
    assert out.columns[-1] == "history"
    row = out.row(0, named=True)
    assert (row["name"], row["version"], row["start"], row["stop"]) == (
        "g1.2",
        2,
        100,
        260,
    )  # URS3 retired: membership changed
