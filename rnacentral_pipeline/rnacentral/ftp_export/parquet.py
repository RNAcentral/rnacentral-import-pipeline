# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import re
from pathlib import Path

import polars as pl

from rnacentral_pipeline.rnacentral.ftp_export.sequences_json import parse

# active.sql is shared with the FTP export, where it is a `COPY (...) TO STDOUT`
# statement. polars' read_database needs a row-returning SELECT, so we unwrap the
# COPY here rather than forking the (shared) .sql file.
_COPY_RE = re.compile(r"COPY\s*\((.*)\)\s*TO\s+STDOUT", re.IGNORECASE | re.DOTALL)


def unwrap_copy(query_str: str) -> str:
    """Return the inner SELECT of a `COPY (...) TO STDOUT`, else the query unchanged."""
    match = _COPY_RE.search(query_str)
    return match.group(1).strip() if match else query_str.strip()


def query_2_dataframe(query_path: str | Path, conn) -> pl.DataFrame:
    """
    Use polars to read a query direct from the database to a dataframe
    expects a psycopg2 connection object
    """
    query_str = Path(query_path).read_text()
    return pl.read_database(unwrap_copy(query_str), conn)


def ndjson_2_dataframe(json_path: str | Path) -> pl.DataFrame:
    """
    Build a dataframe from the active-sequence dump file (the same JSON-lines the
    FASTA export consumes), decoded via the shared parser so the two exports agree.
    """
    # ponytail: materialises the whole dump in memory, matching the DB path; the
    # pipeline runs this on a high-memory node. Switch to a batched pl.read_ndjson
    # (with backslash un-doubling) if memory becomes the limiter.
    with open(json_path) as handle:
        return pl.DataFrame(list(parse(handle)))


if __name__ == "__main__":
    # ponytail: unwrap is the only non-trivial logic here; DB read is exercised in the release run.
    assert unwrap_copy("COPY (\nSELECT 1\n) TO STDOUT").strip() == "SELECT 1"
    assert unwrap_copy("select a from b") == "select a from b"
    print("ok")
