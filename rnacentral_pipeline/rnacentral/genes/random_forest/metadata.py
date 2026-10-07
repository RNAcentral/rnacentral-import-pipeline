# -*- coding: utf-8 -*-

"""
Choosing a description and an RNA type for each gene from the evidence on its
member transcripts.

The evidence is a table with one row per (URS, source record): the accessions
of the URS from the expert databases, plus its Rfam and R2DT hits. Each row
gets a description score and a type score; a gene takes the description of its
highest-scoring row and the type of its highest-scoring row, over all members.

Everything here works on whole tables, with no database access, so a gene costs
nothing beyond the join. Reducing to the best row per URS first is exact (a max
of maxes) and keeps the table small enough to build per batch of URS.

Ties are broken deterministically: rows are ordered by (database, description,
rna_type) first and the first maximum wins.
"""

import typing as ty

import polars as pl

EVIDENCE_SCHEMA = {
    "urs_taxid": pl.Utf8,
    "database": pl.Utf8,
    "description": pl.Utf8,
    "rna_type": pl.Utf8,
    "cm_overlap": pl.Float64,
}

BEST_SCHEMA = {
    "urs_taxid": pl.Utf8,
    "description": pl.Utf8,
    "desc_score": pl.Float64,
    "rna_type": pl.Utf8,
    "type_score": pl.Float64,
}


def _mapped(values: pl.Series, function: ty.Callable, dtype) -> pl.Series:
    """Apply `function` once per distinct value, then fan the result back out."""
    distinct = values.drop_nulls().unique().to_list()
    lookup = pl.DataFrame(
        {"value": distinct, "result": [function(v) for v in distinct]},
        schema={"value": values.dtype, "result": dtype},
    )
    return values.to_frame("value").join(lookup, on="value", how="left")["result"]


def score_evidence(
    evidence: pl.DataFrame,
    description_scores: ty.Mapping[str, float],
    type_scores: ty.Mapping[str, float],
    entropy: ty.Callable[[str], float],
    specificity: ty.Callable[[str], float],
) -> pl.DataFrame:
    """Add desc_score and type_score columns to the evidence table."""
    evidence = evidence.with_columns(pl.col("cm_overlap").fill_null(0.0))
    db_desc = pl.col("database").replace_strict(
        description_scores, default=0, return_dtype=pl.Float64
    )
    db_type = pl.col("database").replace_strict(
        type_scores, default=0, return_dtype=pl.Float64
    )
    return evidence.with_columns(
        desc_score=db_desc
        + _mapped(evidence["description"], entropy, pl.Float64).fill_null(0.0),
        type_score=db_type
        + _mapped(evidence["rna_type"], specificity, pl.Float64).fill_null(-1000.0)
        + pl.col("cm_overlap") * 100.0,
    )


def _best(frame: pl.DataFrame, by: str, keys: ty.List[str]) -> pl.DataFrame:
    return frame.group_by(keys, maintain_order=True).agg(
        pl.col("description").get(pl.col("desc_score").arg_max()).alias("description")
        if by == "desc_score"
        else pl.col("rna_type").get(pl.col("type_score").arg_max()).alias("rna_type"),
        pl.col(by).max(),
    )


def best_per_urs(scored: pl.DataFrame) -> pl.DataFrame:
    """The best description and the best type for each URS, with their scores."""
    scored = scored.sort(
        ["urs_taxid", "database", "description", "rna_type"], nulls_last=True
    )
    desc = _best(scored, "desc_score", ["urs_taxid"])
    rna_type = _best(scored, "type_score", ["urs_taxid"])
    return desc.join(rna_type, on="urs_taxid", how="full", coalesce=True).select(
        list(BEST_SCHEMA)
    )


def choose_gene_metadata(
    per_urs: pl.DataFrame,
    members: pl.DataFrame,
    normalize_type: ty.Callable[[ty.Optional[str]], str],
) -> pl.DataFrame:
    """
    `members` has columns name and urs_taxid (one row per member). Returns
    name, description, rna_type; genes none of whose members have any evidence
    are absent from the result.
    """
    joined = members.join(per_urs, on="urs_taxid", how="inner").sort(
        ["name", "urs_taxid"]
    )
    chosen = joined.group_by("name", maintain_order=True).agg(
        pl.col("description").get(pl.col("desc_score").arg_max()),
        pl.col("rna_type").get(pl.col("type_score").arg_max()),
    )
    # A gene whose rows all have a missing type gets the normaliser's answer for None
    fallback = pl.lit(normalize_type(None), dtype=pl.Utf8)
    return chosen.with_columns(
        rna_type=pl.coalesce(
            _mapped(chosen["rna_type"], normalize_type, pl.Utf8), fallback
        )
    )
