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

# Polars replacement for the Rust `precompute metadata group`+`merge`
# binaries. Reads N raw per-type files (one row per hit, typed Parquet
# produced by `rnac precompute extract-query` - see
# docs/superpowers/specs/2026-09-29-duckdb-parquet-metadata-extraction-design.md)
# and lazily joins them onto a `basic`-derived spine, once. `build()` then
# filters and sinks that one lazy plan once per id range in urs_taxid.csv,
# producing N range-scoped files instead of Rust's one monolithic
# metadata.json.
# See docs/superpowers/specs/2026-09-29-polars-precompute-metadata-build-design.md.

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import polars as pl

LOGGER = logging.getLogger(__name__)


def as_list(raw: pl.LazyFrame, name: str) -> pl.LazyFrame:
    """AnyNumber cardinality: group raw per-hit rows into a list of structs per id."""
    return raw.group_by("id").agg(pl.struct(pl.exclude("id")).alias(name))


def as_optional_struct(raw: pl.LazyFrame, name: str) -> pl.LazyFrame:
    """ZeroOrOne cardinality: collapse at-most-one-row-per-id raw rows into a nullable struct."""
    return raw.group_by("id").agg(pl.struct(pl.exclude("id")).first().alias(name))


def as_optional_scalar(raw: pl.LazyFrame, source_col: str, name: str) -> pl.LazyFrame:
    """ZeroOrOne cardinality: collapse at-most-one-row-per-id raw rows into a nullable scalar."""
    return raw.group_by("id").agg(pl.col(source_col).first().alias(name))


@dataclass
class Attachment:
    filename: str
    join: Callable[[pl.LazyFrame], pl.LazyFrame]
    # "list" (AnyNumber - multiple raw rows per id is normal) or "optional"
    # (ZeroOrOne - Rust's grouper.rs errors loudly on >1 row; this design
    # silently keeps the first via .first() instead, so >1 is logged as a
    # data-quality warning rather than either crashing or vanishing).
    cardinality: str = "list"
    fill_nulls: dict = field(default_factory=dict)


def read_spine(raw_dir: Path) -> pl.LazyFrame:
    return pl.scan_parquet(raw_dir / "basic.parquet")


def build_plan(raw_dir: Path, attachments: list) -> pl.LazyFrame:
    plan = read_spine(raw_dir)
    for attachment in attachments:
        path = raw_dir / attachment.filename
        if not path.exists():
            raise FileNotFoundError(f"missing attachment file: {path}")
        joined = attachment.join(pl.scan_parquet(path))
        plan = plan.join(joined, on="id", how="left")
        if attachment.fill_nulls:
            plan = plan.with_columns(
                pl.col(name).fill_null(default)
                for name, default in attachment.fill_nulls.items()
            )
    return plan


def orf_join(raw: pl.LazyFrame) -> pl.LazyFrame:
    """Port of orf.rs: orf_info.sources = unique non-null `source` values
    across the group (None if empty); possible_orf = the FIRST non-null
    `is_protein_coding` in the group (Rust's find_map, not any())."""
    grouped = raw.group_by("id").agg(
        pl.col("source").drop_nulls().unique().alias("sources"),
        pl.col("is_protein_coding").drop_nulls().first().alias("possible_orf"),
    )
    return grouped.with_columns(
        pl.when(pl.col("sources").list.len() > 0)
        .then(pl.struct(pl.col("sources")))
        .otherwise(None)
        .alias("orf_info")
    ).drop("sources")


ATTACHMENTS = [
    Attachment(
        "coordinates.parquet",
        lambda raw: as_list(raw, "coordinates"),
        fill_nulls={"coordinates": []},
    ),
    Attachment(
        "rfam-hits.parquet",
        lambda raw: as_list(raw, "rfam_hits"),
        fill_nulls={"rfam_hits": []},
    ),
    Attachment(
        "r2dt-hits.parquet",
        lambda raw: as_optional_struct(raw, "r2dt_hits"),
        cardinality="optional",
    ),
    Attachment(
        "previous.parquet",
        lambda raw: as_optional_struct(raw, "previous"),
        cardinality="optional",
    ),
    Attachment("orfs.parquet", orf_join),
    Attachment(
        "stopfree.parquet",
        lambda raw: as_optional_scalar(
            raw, "is_protein_coding", "possible_orf_stopfree"
        ),
        cardinality="optional",
    ),
    Attachment(
        "tcode.parquet",
        lambda raw: as_optional_scalar(raw, "is_protein_coding", "possible_orf_tcode"),
        cardinality="optional",
    ),
]


def read_ranges(path: Path):
    with open(path) as f:
        return [(int(a), int(b), int(c)) for a, b, c in csv.reader(f) if a]


def _log_data_quality_warnings(raw_dir: Path, attachments) -> None:
    # Two things build_plan's left join can't surface on its own, logged
    # here instead of vanishing silently: (1) an attachment id absent from
    # basic - dropped by the join (Rust's dense-fill loop would error
    # loudly instead); (2) for a "optional" (ZeroOrOne) attachment, an id
    # with more than one raw row - Rust's grouper.rs errors loudly
    # ("Too many items... expected 0 or 1"), this design's .first() keeps
    # the first row silently. Cheap: only scans id (+ group sizes), not the
    # full nested multi-attachment plan.
    spine_ids = read_spine(raw_dir).select("id").unique()
    for attachment in attachments:
        raw = pl.scan_parquet(raw_dir / attachment.filename)
        raw_ids = raw.select("id").unique()
        unmatched = (
            raw_ids.join(spine_ids, on="id", how="anti")
            .select(pl.len())
            .collect()
            .item()
        )
        if unmatched:
            LOGGER.warning(
                "%s: %d id(s) not present in basic, dropped",
                attachment.filename,
                unmatched,
            )
        if attachment.cardinality == "optional":
            duplicated = (
                raw.group_by("id")
                .len()
                .filter(pl.col("len") > 1)
                .select(pl.len())
                .collect()
                .item()
            )
            if duplicated:
                LOGGER.warning(
                    "%s: %d id(s) had more than one row for a single-value "
                    "attachment, kept the first",
                    attachment.filename,
                    duplicated,
                )


def build(ranges_path: Path, raw_dir: Path, output_dir: Path) -> None:
    _log_data_quality_warnings(raw_dir, ATTACHMENTS)
    plan = build_plan(raw_dir, ATTACHMENTS)
    ranges = read_ranges(ranges_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = output_dir / "manifest.csv"
    with open(manifest_path, "w", newline="") as manifest_file:
        writer = csv.writer(manifest_file)
        for upi_min, ut_min, ut_max in ranges:
            chunk_path = output_dir / f"metadata-{upi_min}.json"
            plan.filter(pl.col("id").is_between(ut_min, ut_max)).sink_ndjson(chunk_path)
            writer.writerow([upi_min, chunk_path.resolve()])

    LOGGER.info("wrote %d metadata chunks to %s", len(ranges), output_dir)
