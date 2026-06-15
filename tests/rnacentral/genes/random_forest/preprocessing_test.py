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
Functional-equivalence tests for the gene random-forest preprocessing.

The vectorised polars implementation (``polars_work_function`` and friends) is
checked against an independent, deliberately naive row-by-row reference oracle
(``reference_features``) built from the small per-pair helpers. The oracle
encodes the *intended* semantics of the feature generation, so these tests pin
the contract and let the vectorised query be restructured (split up, spilled to
disk, etc.) with confidence that the output is unchanged.
"""

import numpy as np
import polars as pl
import pytest

gpp = pytest.importorskip(
    "gene_preprocessing",
    reason="Rust gene_preprocessing extension not built "
    "(cd utils/genes-preprocessing && maturin develop)",
)

from rnacentral_pipeline.rnacentral.genes.random_forest import preprocessing as P
from rnacentral_pipeline.rnacentral.genes.random_forest.preprocessing import (
    exon_overlap_tup,
)

# A tiny, fully controlled SO "model": unit vectors so dot-product similarities
# are exactly known. A and C are identical (sim 1.0); B is orthogonal (sim 0.0);
# D is at 45 degrees to A/C (sim ~0.707).
SO_MODEL = {
    "SO:A": np.array([1.0, 0.0], dtype=np.float32),
    "SO:B": np.array([0.0, 1.0], dtype=np.float32),
    "SO:C": np.array([1.0, 0.0], dtype=np.float32),
    "SO:D": np.array([1.0, 1.0], dtype=np.float32) / np.sqrt(2.0),
}


@pytest.fixture(autouse=True)
def _inject_so_model(monkeypatch):
    """Make get_type_similarity use the controlled model, no Word2Vec/IO."""
    monkeypatch.setattr(P, "_SO_MODEL", SO_MODEL)
    # get_type_similarity is lru_cached; clear so values don't leak across tests.
    P.get_type_similarity.cache_clear()
    yield
    P.get_type_similarity.cache_clear()


def make_transcripts(rows):
    """Build a transcripts DataFrame matching the production schema."""
    return pl.DataFrame(
        rows,
        orient="row",
        schema={
            "region_name": pl.Utf8,
            "region_start": pl.Int64,
            "region_stop": pl.Int64,
            "chromosome": pl.Utf8,
            "assembly_id": pl.Utf8,
            "strand": pl.Int64,
            "so_type": pl.Utf8,
            "exon_start": pl.List(pl.Int64),
            "exon_stop": pl.List(pl.Int64),
        },
    )


def reference_features(transcripts, nearby_distance=1000):
    """
    Independent row-by-row reference for the feature table.

    Mirrors the join predicates (same chromosome / assembly / strand, within
    nearby_distance, lower region_start is the 'a' member) and computes each
    feature with the simple per-pair helpers rather than the vectorised query.
    """
    recs = transcripts.sort("region_start").to_dicts()
    out = []
    for i, a in enumerate(recs):
        for j in range(i + 1, len(recs)):
            b = recs[j]
            if a["chromosome"] != b["chromosome"]:
                continue
            if a["assembly_id"] != b["assembly_id"]:
                continue
            if a["strand"] != b["strand"]:
                continue
            if abs(b["region_start"] - a["region_start"]) > nearby_distance:
                continue

            if a["strand"] == 1:
                ex5a = (a["exon_start"][0], a["exon_stop"][0])
                ex5b = (b["exon_start"][0], b["exon_stop"][0])
            else:
                ex5a = (a["exon_start"][-1], a["exon_stop"][-1])
                ex5b = (b["exon_start"][-1], b["exon_stop"][-1])

            overlap = exon_overlap_tup(ex5a, ex5b)
            dta, dta_3p = gpp.distance_to_agreement(list(ex5a), list(ex5b))

            exons_a = list(zip(a["exon_start"], a["exon_stop"]))
            exons_b = list(zip(b["exon_start"], b["exon_stop"]))
            count = gpp.count_overlap(exons_a, exons_b, 0.9)

            type_sim = float(
                np.dot(SO_MODEL[a["so_type"]], SO_MODEL[b["so_type"]])
            )

            out.append(
                {
                    "comparison": f"{a['region_name']} vs {b['region_name']}",
                    "5p_exon_overlap": float(overlap),
                    "5p_exon_dta": int(dta),
                    "5p_exon_3p_dta": int(dta_3p),
                    "exons_overlapping": int(count),
                    "strand": 1,  # join keeps only same-strand pairs
                    "type_sim": type_sim,
                }
            )
    return out


def as_lookup(records):
    """Index feature rows by their comparison string for order-free compare."""
    return {r["comparison"]: r for r in records}


def assert_features_match(produced_df, expected_records):
    """Compare a produced feature frame against reference records, order-free."""
    produced = as_lookup(produced_df.to_dicts())
    expected = as_lookup(expected_records)

    assert set(produced) == set(expected), "different set of comparisons"

    for key, exp in expected.items():
        got = produced[key]
        assert got["5p_exon_dta"] == exp["5p_exon_dta"], key
        assert got["5p_exon_3p_dta"] == exp["5p_exon_3p_dta"], key
        assert got["exons_overlapping"] == exp["exons_overlapping"], key
        assert got["strand"] == exp["strand"], key
        assert got["5p_exon_overlap"] == pytest.approx(exp["5p_exon_overlap"]), key
        assert got["type_sim"] == pytest.approx(exp["type_sim"], abs=1e-6), key


@pytest.fixture
def rich_chromosome():
    """
    A single chromosome exercising: a partial 5' overlap, a full exon overlap,
    a non-overlapping near pair, an out-of-range pair, an opposite-strand pair,
    and varied type similarities.
    """
    return make_transcripts(
        [
            # name, start, stop, chrom, asm, strand, so_type, exon_start, exon_stop
            ("a", 100, 400, "1", "GRCh38", 1, "SO:A", [100, 300], [200, 400]),
            ("b", 150, 450, "1", "GRCh38", 1, "SO:C", [150, 350], [250, 450]),
            ("c", 700, 900, "1", "GRCh38", 1, "SO:D", [700], [900]),
            ("d", 5000, 5200, "1", "GRCh38", 1, "SO:A", [5000], [5200]),
            ("e", 160, 460, "1", "GRCh38", -1, "SO:B", [160, 360], [260, 460]),
        ]
    )


def test_work_function_matches_reference(rich_chromosome):
    produced = P.polars_work_function(rich_chromosome, nearby_distance=1000)
    expected = reference_features(rich_chromosome, nearby_distance=1000)
    assert_features_match(produced, expected)


def test_work_function_output_schema(rich_chromosome):
    produced = P.polars_work_function(rich_chromosome, nearby_distance=1000)
    # Same columns and dtypes as the canonical empty feature frame.
    assert produced.schema == P.empty_features.schema


def test_opposite_strand_pairs_excluded(rich_chromosome):
    produced = P.polars_work_function(rich_chromosome, nearby_distance=1000)
    comparisons = set(produced["comparison"].to_list())
    # 'e' is on the - strand and near a/b but must never pair with them.
    assert not any("e" in c for c in comparisons)


def test_distance_filtering(rich_chromosome):
    near = set(
        P.polars_work_function(rich_chromosome, nearby_distance=1000)[
            "comparison"
        ].to_list()
    )
    # 'c' (start 700) is within 1000 of 'a' (start 100) -> paired.
    assert "a vs c" in near
    # 'd' (start 5000) is 4900 away -> excluded at distance 1000.
    assert "a vs d" not in near
    wide = set(
        P.polars_work_function(rich_chromosome, nearby_distance=10000)[
            "comparison"
        ].to_list()
    )
    # ...but included once the threshold is wide enough (4900 < 10000).
    assert "a vs d" in wide


def test_exon_overlap_threshold():
    # x and y have two near-identical exons (>90% overlap) -> both counted.
    # z's exons are offset by half their width (50% overlap) -> not counted.
    transcripts = make_transcripts(
        [
            ("x", 100, 400, "1", "GRCh38", 1, "SO:A", [100, 300], [200, 400]),
            ("y", 105, 405, "1", "GRCh38", 1, "SO:A", [105, 305], [205, 405]),
            ("z", 150, 450, "1", "GRCh38", 1, "SO:A", [150, 350], [250, 450]),
        ]
    )
    produced = as_lookup(
        P.polars_work_function(transcripts, nearby_distance=1000).to_dicts()
    )
    assert produced["x vs y"]["exons_overlapping"] == 2
    assert produced["x vs z"]["exons_overlapping"] == 0
    # Cross-check against the independent oracle.
    assert_features_match(
        P.polars_work_function(transcripts, nearby_distance=1000),
        reference_features(transcripts, nearby_distance=1000),
    )


def test_polars_preprocessing_matches_reference_multichrom():
    transcripts = make_transcripts(
        [
            ("a", 100, 400, "1", "GRCh38", 1, "SO:A", [100, 300], [200, 400]),
            ("b", 150, 450, "1", "GRCh38", 1, "SO:C", [150, 350], [250, 450]),
            ("c", 100, 300, "2", "GRCh38", 1, "SO:D", [100], [300]),
            ("d", 200, 400, "2", "GRCh38", 1, "SO:A", [200], [400]),
            ("lonely", 999, 1100, "3", "GRCh38", 1, "SO:A", [999], [1100]),
        ]
    )
    produced = P.polars_preprocessing(transcripts, nearby_distance=1000)
    expected = reference_features(transcripts, nearby_distance=1000)
    assert_features_match(produced, expected)
    # chromosome 3 has a single transcript -> contributes no pairs.
    assert all("lonely" not in c for c in produced["comparison"].to_list())


def test_polars_preprocessing_label_column():
    transcripts = make_transcripts(
        [
            ("a", 100, 400, "1", "GRCh38", 1, "SO:A", [100], [400]),
            ("b", 150, 450, "1", "GRCh38", 1, "SO:C", [150], [450]),
        ]
    )
    produced = P.polars_preprocessing(transcripts, nearby_distance=1000, label=1)
    assert "label" in produced.columns
    assert produced["label"].to_list() == [1]


def test_sink_preprocessing_matches_in_memory(tmp_path):
    transcripts = make_transcripts(
        [
            ("a", 100, 400, "1", "GRCh38", 1, "SO:A", [100, 300], [200, 400]),
            ("b", 150, 450, "1", "GRCh38", 1, "SO:C", [150, 350], [250, 450]),
            ("c", 100, 300, "2", "GRCh38", 1, "SO:D", [100], [300]),
            ("d", 200, 400, "2", "GRCh38", 1, "SO:A", [200], [400]),
        ]
    )
    in_memory = P.polars_preprocessing(transcripts, nearby_distance=1000)

    out = tmp_path / "nested" / "features.parquet"
    n = P.sink_preprocessing(transcripts, out, nearby_distance=1000)

    assert out.exists()
    assert n == in_memory.height
    sunk = pl.read_parquet(out)
    assert sunk.schema == in_memory.schema
    assert_features_match(sunk, reference_features(transcripts, nearby_distance=1000))


def test_empty_input_returns_empty_schema():
    transcripts = make_transcripts([])
    produced = P.polars_preprocessing(transcripts, nearby_distance=1000)
    assert produced.height == 0
    assert produced.schema == P.empty_features.schema


def test_sink_empty_input(tmp_path):
    transcripts = make_transcripts([])
    out = tmp_path / "features.parquet"
    n = P.sink_preprocessing(transcripts, out, nearby_distance=1000)
    assert n == 0
    assert pl.read_parquet(out).height == 0
