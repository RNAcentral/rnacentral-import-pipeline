# -*- coding: utf-8 -*-

import polars as pl

from rnacentral_pipeline.rnacentral.genes.random_forest import metadata

DESC = {"hgnc": 100, "ena": 2}
TYPE = {"rfam": 100, "ena": 0}
SPECIFICITY = {"SO:miRNA": 5.0, "SO:ncRNA": 1.0}
ROOT = "SO:ncRNA"


def evidence(rows):
    return pl.DataFrame(rows, schema=metadata.EVIDENCE_SCHEMA, orient="row")


def score(frame):
    return metadata.score_evidence(
        frame,
        DESC,
        TYPE,
        entropy=lambda d: float(len(set(d))),
        specificity=SPECIFICITY.get,
    )


def choose(rows, genes):
    per_urs = metadata.best_per_urs(score(evidence(rows)))
    members = pl.DataFrame(
        genes, schema={"name": pl.Utf8, "urs_taxid": pl.Utf8}, orient="row"
    )
    out = metadata.choose_gene_metadata(
        per_urs, members, lambda t: t if t in SPECIFICITY else ROOT
    )
    return {r["name"]: (r["description"], r["rna_type"]) for r in out.to_dicts()}


def test_description_and_type_come_from_the_best_scoring_rows_across_members():
    rows = [
        ("U1_9", "ena", "vague", "SO:ncRNA", 0.0),
        ("U2_9", "hgnc", "gene symbol", "SO:miRNA", 0.0),
        ("U2_9", "rfam", "family model", "SO:ncRNA", 1.0),
    ]
    got = choose(rows, [("g1", "U1_9"), ("g1", "U2_9")])
    # hgnc (+100) wins the description; rfam (100 + specificity 1 + overlap 100) beats ena and miRNA-by-hgnc (0 + 5)
    assert got == {"g1": ("gene symbol", "SO:ncRNA")}


def test_best_per_urs_is_the_same_as_scoring_all_rows_at_once():
    rows = [
        ("U1_9", "ena", "a", "SO:ncRNA", 0.0),
        ("U1_9", "hgnc", "bc", "SO:miRNA", 0.0),
        ("U2_9", "ena", "d", "SO:ncRNA", 0.0),
    ]
    whole = choose(rows, [("g", "U1_9"), ("g", "U2_9")])
    parts = pl.concat(
        [
            metadata.best_per_urs(score(evidence(rows[:2]))),
            metadata.best_per_urs(score(evidence(rows[2:]))),
        ]
    )
    members = pl.DataFrame(
        [("g", "U1_9"), ("g", "U2_9")],
        schema={"name": pl.Utf8, "urs_taxid": pl.Utf8},
        orient="row",
    )
    batched = metadata.choose_gene_metadata(parts, members, lambda t: t or ROOT)
    assert (batched["description"][0], batched["rna_type"][0]) == whole["g"]


def test_genes_without_evidence_are_absent_and_ties_are_deterministic():
    rows = [
        ("U1_9", "ena", "same", "SO:ncRNA", 0.0),
        ("U1_9", "ena", "same", "SO:ncRNA", 0.0),
    ]
    assert choose(rows, [("g1", "U1_9"), ("g2", "U9_9")]) == {
        "g1": ("same", "SO:ncRNA")
    }


def test_unknown_types_and_null_overlap_are_handled():
    rows = [("U1_9", "ena", "x", None, None), ("U1_9", "rfam", "y", "SO:odd", None)]
    got = choose(rows, [("g", "U1_9")])
    assert got["g"][1] == ROOT  # unknown type normalised to the generic term
