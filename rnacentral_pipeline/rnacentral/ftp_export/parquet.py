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

from itertools import islice
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq

from rnacentral_pipeline.rnacentral.ftp_export.sequences_json import parse


def _records_2_parquet(records, output: str | Path, batch_size: int):
    """
    Write an iterator of dicts to parquet a batch at a time, so peak memory is one
    batch. active.sql has a fixed schema (a single json_build_object), so the first
    batch defines the schema and later batches cast to match.
    """
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


def ndjson_2_parquet(
    json_path: str | Path, output: str | Path, batch_size: int = 50_000
):
    """
    Stream the active-sequence dump (the same JSON-lines the FASTA export consumes,
    decoded via the shared parser) into parquet a batch at a time.
    """
    with open(json_path) as handle:
        _records_2_parquet(parse(handle), output, batch_size)


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        # streaming write across multiple batches must round-trip every row
        src = Path(tmp) / "dump.json"
        out = Path(tmp) / "out.parquet"
        src.write_text(
            "".join(f'{{"id": "u{i}", "sequence": "ACGT"}}\n' for i in range(5))
        )
        ndjson_2_parquet(src, out, batch_size=2)
        assert pl.read_parquet(out).height == 5, pl.read_parquet(out)

    print("ok")
