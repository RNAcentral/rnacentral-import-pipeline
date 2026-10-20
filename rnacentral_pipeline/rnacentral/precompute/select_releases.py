# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

# Polars port of utils/precompute/src/releases.rs::select_new - selects
# which urs_taxid pairs need recomputing, by comparing each pair's xref
# release against its last-known-precomputed release. Replaces Rust's
# sorted-merge outer join (which needed both inputs pre-sorted, hence the
# external `sort` calls in workflows/precompute/build_urs_table.nf) with a
# plain Polars hash join, which doesn't.
#
# xref_release is null only if urs_taxid isn't in xref at all (impossible
# by construction - it comes from xref); known_release is null if the pair
# has never been precomputed. See
# docs/superpowers/plans/2026-09-29-precompute-select-migration.md for the
# full design and the Rust test cases this was verified against.
#
# Lazy throughout (scan_parquet + collect(engine="streaming") only at the
# end), not eager read_parquet - measured directly at a representative
# scale (20M/18M synthetic rows): eager peaked at ~4.96GB RSS, lazy+streaming
# at ~3.45GB and ran ~2.4x faster (projection pushdown means only
# urs_taxid, not both release columns, is ever fully materialized). xref
# is expected to be real-DB-scale (tens of millions of rows), so this
# isn't a style preference.

from pathlib import Path

import polars as pl


def select_new(xref_lf: pl.LazyFrame, known_lf: pl.LazyFrame) -> list[str]:
    joined = xref_lf.join(known_lf, on="urs_taxid", how="full", coalesce=True)

    too_new = joined.filter(
        pl.col("known_release").is_not_null()
        & pl.col("xref_release").is_not_null()
        & (pl.col("known_release") > pl.col("xref_release"))
    ).collect(engine="streaming")
    if too_new.height > 0:
        raise ValueError(
            f"{too_new.height} urs_taxid have a known/precompute release "
            "newer than xref (should never happen): "
            f"{too_new['urs_taxid'].to_list()[:10]}"
        )

    selected = (
        joined.filter(
            pl.col("xref_release").is_not_null()
            & (
                pl.col("known_release").is_null()
                | (pl.col("xref_release") > pl.col("known_release"))
            )
        )
        .select("urs_taxid")
        .collect(engine="streaming")
    )
    return selected["urs_taxid"].to_list()


def write(xref_path: Path, known_path: Path, output_path: Path) -> None:
    xref_lf = pl.scan_parquet(xref_path)
    known_lf = pl.scan_parquet(known_path)
    selected = select_new(xref_lf, known_lf)
    with open(output_path, "w") as f:
        for urs_taxid in selected:
            f.write(urs_taxid + "\n")
