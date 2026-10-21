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

# Runs one of files/precompute/queries/*.sql (or get-accessions/query.sql)
# against Postgres via Polars' connectorx engine and writes the typed
# result straight to Parquet - replaces psql piping json_build_object(...)
# output to a file. See
# docs/superpowers/specs/2026-09-29-duckdb-parquet-metadata-extraction-design.md.
#
# :min/:max are substituted as literal text before the query reaches
# connectorx (which has no named-parameter support via
# pl.read_database_uri) - the same mechanism psql's `-v min=$min` already
# uses today, not a new kind of risk. Every rewritten query file has its
# own `WHERE (:min IS NULL OR <id_column> BETWEEN :min AND :max)` clause,
# so an unranged call just substitutes NULL for both.

from pathlib import Path

import polars as pl


def _substitute_range(sql: str, range: tuple[int, int] | None) -> str:
    if range is None:
        min_text = max_text = "NULL"
    else:
        min_text, max_text = str(int(range[0])), str(int(range[1]))
    return sql.replace(":min", min_text).replace(":max", max_text)


def extract_query(
    sql_path: Path,
    output_path: Path,
    pg_uri: str,
    range: tuple[int, int] | None = None,
) -> None:
    sql = sql_path.read_text()
    query = _substitute_range(sql, range)
    df = pl.read_database_uri(query, pg_uri, engine="connectorx")
    df.write_parquet(output_path)
