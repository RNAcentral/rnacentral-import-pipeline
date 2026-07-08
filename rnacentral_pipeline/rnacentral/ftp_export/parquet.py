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

from rnacentral_pipeline.rnacentral.ftp_export.sequences_json import decode_entry, parse


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


def copy_query_2_parquet(
    query_path: str | Path, conn, output: str | Path, batch_size: int = 250_000
):
    """
    Stream a `COPY (...) TO STDOUT` query straight into parquet. Unlike read_database
    this never materialises the whole (multi-hundred-GB) result set in memory, and
    unlike the FASTA path it needs no intermediate dump file: psycopg2 pushes chunks
    to sink.write(), which splits them into JSON lines, decodes them through the same
    shared parser, and flushes a batch at a time. Only `output` is ever written.
    """
    query_str = Path(query_path).read_text()
    sink = _ParquetLineSink(output, batch_size)
    with conn.cursor() as cur:
        cur.copy_expert(query_str, sink)
    sink.close()


class _ParquetLineSink:
    """File-like sink for psycopg2 copy_expert: buffers partial lines across chunks
    and drives the batched ParquetWriter. Mirrors _records_2_parquet's writer loop
    but push-driven (copy_expert calls write()), so a threaded pipe isn't needed."""

    def __init__(self, output: str | Path, batch_size: int):
        self.output = str(output)
        self.batch_size = batch_size
        self.buf = ""
        self.batch: list[dict] = []
        self.writer = None

    def write(self, data):
        if isinstance(data, bytes):
            data = data.decode()
        self.buf += data
        while (nl := self.buf.find("\n")) != -1:
            line, self.buf = self.buf[:nl], self.buf[nl + 1 :]
            if line:
                self.batch.append(decode_entry(line))
                if len(self.batch) >= self.batch_size:
                    self._flush()

    def _flush(self):
        if not self.batch:
            return
        table = pl.DataFrame(self.batch).to_arrow()
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output, table.schema)
        self.writer.write_table(table.cast(self.writer.schema))
        self.batch = []

    def close(self):
        if self.buf.strip():
            self.batch.append(decode_entry(self.buf.strip()))
        self._flush()
        if self.writer is not None:
            self.writer.close()


def ndjson_2_parquet(json_path: str | Path, output: str | Path, batch_size: int = 50_000):
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
        src.write_text("".join(f'{{"id": "u{i}", "sequence": "ACGT"}}\n' for i in range(5)))
        ndjson_2_parquet(src, out, batch_size=2)
        assert pl.read_parquet(out).height == 5, pl.read_parquet(out)

        # the copy sink must reassemble lines split across write() chunks
        sink_out = Path(tmp) / "sink.parquet"
        sink = _ParquetLineSink(sink_out, batch_size=2)
        blob = "".join(f'{{"id": "u{i}", "sequence": "ACGT"}}\n' for i in range(5))
        for i in range(0, len(blob), 7):  # arbitrary chunking, splits mid-line
            sink.write(blob[i : i + 7])
        sink.close()
        assert pl.read_parquet(sink_out).height == 5, pl.read_parquet(sink_out)
    print("ok")
