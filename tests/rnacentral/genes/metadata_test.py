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

import networkx as nx
import polars as pl
import pytest

from rnacentral_pipeline.rnacentral.genes import merge
from rnacentral_pipeline.rnacentral.genes.metadata import (
    EVIDENCE_SCHEMA,
    compute_metadata,
)

MISC, MIRNA, PRE = "SO:0000673", "SO:0000276", "SO:0001244"


@pytest.fixture(autouse=True)
def tiny_so_graph(monkeypatch):
    g = nx.DiGraph()
    g.add_edges_from([(MIRNA, MISC), (PRE, MIRNA), ("SO:0000655", MISC)])
    monkeypatch.setattr(merge, "get_so_graph", lambda: g)


def members(**genes):
    return pl.DataFrame(
        [
            {"internal_name": g, "region_name": f"{u}@1/1-2:+"}
            for g, urs in genes.items()
            for u in urs
        ]
    )


def ev(*rows):
    return pl.DataFrame(rows, schema=EVIDENCE_SCHEMA, orient="row")


def test_best_description_and_type_per_gene():
    m, missing = compute_metadata(
        members(g1=["URS1_9606"]),
        ev(
            ("URS1_9606", "ENA", "Homo sapiens RNA", "SO:0000655", 0.0),
            (
                "URS1_9606",
                "MIRBASE",
                "Homo sapiens (human) microRNA hsa-mir-1",
                MIRNA,
                0.0,
            ),
        ),
    )
    assert missing.is_empty()
    row = m.row(0, named=True)
    assert row["description"].startswith("Homo sapiens (human) microRNA")
    assert row["rna_type"] == MIRNA
    assert "Homo sapiens" not in row["short_description"]


def test_gene_without_evidence_is_reported_missing():
    m, missing = compute_metadata(
        members(g1=["URS1_9606"], g2=["URS2_9606"]),
        ev(("URS1_9606", "ENA", "x", MIRNA, 0.0)),
    )
    assert m["internal_name"].to_list() == ["g1"]
    assert missing["internal_name"].to_list() == ["g2"]


def test_unknown_or_null_so_id_becomes_generic_ncrna():
    m, _ = compute_metadata(
        members(g1=["URS1_9606"]),
        ev(("URS1_9606", "ENA", "thing", None, 0.0)),
    )
    assert m["rna_type"].to_list() == ["SO:0000655"]


def test_cm_overlap_raises_type_score():
    m, _ = compute_metadata(
        members(g1=["URS1_9606"]),
        ev(
            ("URS1_9606", "RFAM", "model a", MIRNA, 0.1),
            ("URS1_9606", "RFAM", "model b", PRE, 0.9),
        ),
    )
    assert m["rna_type"].to_list() == [PRE]


def test_ties_are_deterministic_whatever_the_row_order():
    rows = [
        ("URS1_9606", "RFAM", "same entropy aa", MIRNA, 0.0),
        ("URS1_9606", "RFAM", "same entropy bb", MIRNA, 0.0),
    ]
    a, _ = compute_metadata(members(g1=["URS1_9606"]), ev(*rows))
    b, _ = compute_metadata(members(g1=["URS1_9606"]), ev(*reversed(rows)))
    assert a.equals(b)


def test_evidence_is_pooled_over_all_members():
    m, _ = compute_metadata(
        members(g1=["URS1_9606", "URS2_9606"]),
        ev(
            ("URS1_9606", "ENA", "vague", "SO:0000655", 0.0),
            ("URS2_9606", "HGNC", "Homo sapiens (human) MIR1 microRNA 1", MIRNA, 0.0),
        ),
    )
    assert "MIR1" in m["description"][0]


def test_duplicate_member_urs_does_not_double_count():
    two_regions = pl.DataFrame(
        {
            "internal_name": ["g1", "g1"],
            "region_name": ["URS1_9606@1/1-2:+", "URS1_9606@1/5-6:+"],
        }
    )
    m, _ = compute_metadata(two_regions, ev(("URS1_9606", "ENA", "x", MIRNA, 0.0)))
    assert m.height == 1
