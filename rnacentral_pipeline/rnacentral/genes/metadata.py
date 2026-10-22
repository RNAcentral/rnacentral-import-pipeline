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
Gene description and RNA type, computed without touching the database.

The evidence (accessions, Rfam and R2DT hits for the taxon's URS_taxids) is
extracted once in bulk; `compute_metadata` then scores it in polars and picks the
best description and RNA type per gene.
"""

import polars as pl

from rnacentral_pipeline.rnacentral.genes import merge
from rnacentral_pipeline.rnacentral.precompute import utils

EVIDENCE_SCHEMA = {
    "urs_taxid": pl.Utf8,
    "database": pl.Utf8,
    "description": pl.Utf8,
    "rna_type": pl.Utf8,
    "cm_overlap": pl.Float64,
}

METADATA_SCHEMA = {
    "internal_name": pl.Utf8,
    "description": pl.Utf8,
    "rna_type": pl.Utf8,
    "short_description": pl.Utf8,
}


def _lookup(values, func):
    """Apply func once per distinct value, then map back (entropy and graph walks are slow)."""
    mapping = {v: func(v) for v in values.drop_nulls().unique().to_list()}
    return mapping


def compute_metadata(gene_members, evidence):
    """
    gene_members - frame with internal_name and region_name (one row per member)
    evidence     - frame with EVIDENCE_SCHEMA

    Returns (metadata, missing): metadata has METADATA_SCHEMA, one row per gene that
    has any evidence; missing lists the genes that have none. Ties are broken on the
    description / rna_type text so the result does not depend on row order.
    """
    pairs = gene_members.select(
        "internal_name",
        pl.col("region_name").str.split("@").list.first().alias("urs_taxid"),
    ).unique()
    ev = evidence.cast(EVIDENCE_SCHEMA).with_columns(
        pl.col("cm_overlap").fill_null(0.0), pl.col("description").fill_null("")
    )

    entropy = _lookup(ev["description"], utils.entropy)
    specificity = _lookup(ev["rna_type"], merge.calculate_type_specificity)
    generic = merge.normalize_so_type(None)

    ev = ev.with_columns(
        (
            pl.col("database").replace_strict(merge.description_scores, default=0)
            + pl.col("description").replace_strict(entropy, default=0.0)
        ).alias("desc_score"),
        (
            pl.col("database").replace_strict(merge.type_scores, default=0)
            + pl.col("rna_type").replace_strict(specificity, default=-1000)
            + pl.col("cm_overlap") * 100.0
        ).alias("type_score"),
    )
    joined = pairs.join(ev, on="urs_taxid", how="inner")

    best_desc = (
        joined.sort(
            ["internal_name", "desc_score", "description"],
            descending=[False, True, False],
        )
        .unique(subset="internal_name", keep="first", maintain_order=True)
        .select("internal_name", "description")
    )
    best_type = (
        joined.sort(
            ["internal_name", "type_score", "rna_type"],
            descending=[False, True, False],
            nulls_last=True,
        )
        .unique(subset="internal_name", keep="first", maintain_order=True)
        .select("internal_name", "rna_type")
    )
    metadata = best_desc.join(best_type, on="internal_name")
    ## A missing or unknown SO id collapses to the generic ncRNA term
    normalized = _lookup(metadata["rna_type"], merge.normalize_so_type)
    shorts = _lookup(metadata["description"], merge.remove_species_mention)
    metadata = metadata.with_columns(
        pl.col("rna_type").replace_strict(normalized, default=generic),
        pl.col("description")
        .replace_strict(shorts, default="")
        .alias("short_description"),
    ).sort("internal_name")

    missing = (
        pairs.select("internal_name")
        .unique()
        .join(metadata.select("internal_name"), on="internal_name", how="anti")
        .sort("internal_name")
    )
    return metadata.select(list(METADATA_SCHEMA)), missing
