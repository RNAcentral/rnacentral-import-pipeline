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

# Polars replacement for the Rust `precompute normalize` binary.
#
# Reads the raw, ungrouped accessions file (the shape
# files/precompute/get-accessions/query.sql's psql output has, BEFORE
# Rust's `precompute group-accessions` step - that step is no longer
# needed, this module groups the raw rows itself) and a merged
# metadata.json (still produced by the Rust `precompute metadata merge`
# binary - unchanged for now), and produces the same merged/normalized
# output the Rust `normalize` binary does.
#
# See docs/superpowers/specs/2026-09-28-polars-precompute-normalize-design.md
# for the full design history, including two real bugs found and fixed
# while building this (a Polars schema-inference gotcha on
# sparsely-populated columns, and a psql COPY double-escaping artifact)
# and a documented, deliberate divergence from Rust's actual production
# output: Rust's `normalize` binary applies its backslash-unescape a
# second, redundant time on top of `group-accessions`'s correct pass,
# which corrupts values containing a literal backslash immediately
# followed by a JSON escape-trigger character (t, n, r, b, f, /, \, ").
# This module unescapes exactly once, so it produces the semantically
# correct value there rather than reproducing that corruption.

import json
import logging
import tempfile
from pathlib import Path

import polars as pl

LOGGER = logging.getLogger(__name__)


def _unescape_doubled_backslashes(src: Path, dest: Path) -> None:
    # psql's COPY-to-JSON output double-escapes backslashes: a value that
    # should decode to one literal backslash comes out of the raw dump
    # needing an extra text-level unescape pass before a standard JSON
    # parse gives the right string. Rust's PsqlJsonIterator::next()
    # (rnc-core/psql.rs) does this with buf.replace("\\\\", "\\") on each
    # line before parsing; do the same single pass here, streaming
    # line-by-line so a multi-GB accessions file never gets held whole in
    # memory.
    with open(src) as inp, open(dest, "w") as out:
        for line in inp:
            out.write(line.replace("\\\\", "\\"))


def read_accessions(path: Path) -> pl.DataFrame:
    # Read the raw, ungrouped accession rows (one per xref, the shape
    # get-accessions/query.sql actually produces) and group them ourselves.
    # Reading Rust's pre-grouped {"Multiple": {id, data: [...]}} format
    # instead measured 7GB+ peak memory for a single 25,000-id chunk - the
    # JSON parser over-allocates badly for that doubly-nested
    # List[Struct[List[String]]] shape. Flat rows + group_by measured ~6x
    # less. infer_schema_length=None for the same reason as read_metadata:
    # several fields (e.g. locus_tag) are null for long stretches.
    with tempfile.NamedTemporaryFile(suffix=".json") as tmp:
        unescaped_path = Path(tmp.name)
        _unescape_doubled_backslashes(path, unescaped_path)
        df = pl.read_ndjson(unescaped_path, infer_schema_length=None)
    return df.group_by("id").agg(pl.struct(pl.exclude("id")).alias("data"))


def read_metadata(path: Path) -> pl.DataFrame:
    # infer_schema_length=None scans the whole file for type inference
    # instead of just the first 100 rows (the default). orf_info is null
    # for almost all rows and only occasionally a real struct (CPAT-flagged
    # sequences), so sampling can lock in Null as its type and then crash on
    # the first real value seen later in the file.
    return pl.read_ndjson(path, infer_schema_length=None)


def join_accessions_metadata(
    accessions_df: pl.DataFrame, metadata_df: pl.DataFrame
) -> pl.DataFrame:
    return accessions_df.join(metadata_df, on="id", how="inner")


_ACCESSION_FIELDS = (
    "urs_taxid",
    "accession",
    "is_active",
    "last_release",
    "description",
    "gene",
    "optional_id",
    "database",
    "species",
    "common_name",
    "feature_name",
    "ncrna_class",
    "locus_tag",
    "organelle",
    "lineage",
    "so_rna_type",
)


def _normalize_accession(raw: dict) -> dict:
    result = {field: raw[field] for field in _ACCESSION_FIELDS}
    result["all_species"] = [s for s in raw["all_species"] if s is not None]
    result["all_common_names"] = [s for s in raw["all_common_names"] if s is not None]
    return result


def normalize_row(row: dict) -> dict:
    accessions = [_normalize_accession(a) for a in row["data"]]
    return {
        "urs": row["upi"],
        "taxid": row["taxid"],
        "length": row["length"],
        "last_release": max(a["last_release"] for a in row["data"]),
        "coordinates": row["coordinates"],
        "accessions": accessions,
        "deleted": all(not a["is_active"] for a in row["data"]),
        "previous": row["previous"],
        "rfam_hits": row["rfam_hits"],
        "r2dt_hits": [row["r2dt_hits"]] if row["r2dt_hits"] is not None else [],
        "orf_info": row["orf_info"],
        "possible_orf": row["possible_orf"],
        "possible_orf_stopfree": row["possible_orf_stopfree"],
        "possible_orf_tcode": row["possible_orf_tcode"],
    }


def write_output(rows, path: Path) -> int:
    count = 0
    with open(path, "w") as f:
        for row in rows:
            f.write(json.dumps(row))
            f.write("\n")
            count += 1
    return count


def write(accessions_path: Path, metadata_path: Path, output_path: Path) -> None:
    accessions_df = read_accessions(accessions_path)
    metadata_df = read_metadata(metadata_path)

    joined = join_accessions_metadata(accessions_df, metadata_df)
    accessions_unmatched = accessions_df.height - joined.height
    metadata_unmatched = metadata_df.height - joined.height

    written = write_output(
        (normalize_row(row) for row in joined.iter_rows(named=True)), output_path
    )

    # Grouping raw rows (read_accessions) can never produce an id with an
    # empty accession list - an id with zero accessions just never appears
    # as a group. So unlike the Rust binary's separate "Skipped N urs_taxid
    # with no accessions" count, that case is folded into
    # metadata_unmatched below (a metadata id with no matching accessions
    # group, for whatever reason, including having none at all).
    LOGGER.info(
        "wrote %d rows (accessions unmatched to metadata: %d, metadata "
        "unmatched to accessions, including ids with no accessions at all: %d)",
        written,
        accessions_unmatched,
        metadata_unmatched,
    )
