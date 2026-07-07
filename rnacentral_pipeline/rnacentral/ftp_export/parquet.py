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
from itertools import islice
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq

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


def ndjson_2_parquet(json_path: str | Path, output: str | Path, batch_size: int = 50_000):
    """
    Stream the active-sequence dump (the same JSON-lines the FASTA export consumes,
    decoded via the shared parser) into parquet a batch at a time, so peak memory is
    one batch rather than the whole dump. The dump has a fixed schema (a single
    json_build_object in active.sql), so the first batch defines it and later
    batches cast to match.
    """
    with open(json_path) as handle:
        records = parse(handle)
        writer = None
        try:
            while batch := list(islice(records, batch_size)):
                table = pl.DataFrame(batch).to_arrow()
                if writer is None:
                    writer = pq.ParquetWriter(str(output), table.schema)
                writer.write_table(table.cast(writer.schema))
        finally:
            if writer is not None:
                writer.close()


if __name__ == "__main__":
    import tempfile

    # ponytail: unwrap is the only non-trivial logic here; DB read is exercised in the release run.
    assert unwrap_copy("COPY (\nSELECT 1\n) TO STDOUT").strip() == "SELECT 1"
    assert unwrap_copy("select a from b") == "select a from b"

    # streaming write across multiple batches must round-trip every row
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "dump.json"
        out = Path(tmp) / "out.parquet"
        src.write_text("".join(f'{{"id": "u{i}", "sequence": "ACGT"}}\n' for i in range(5)))
        ndjson_2_parquet(src, out, batch_size=2)
        assert pl.read_parquet(out).height == 5, pl.read_parquet(out)
    print("ok")
