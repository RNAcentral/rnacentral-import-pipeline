# -*- coding: utf-8 -*-

"""
Shared streaming Parquet writer helpers.

Auxiliary parsers and the central :mod:`rnacentral_pipeline.writers` module
both need the same row-buffer-then-flush behaviour: accept positional row
tuples (matching the existing ``csv.writer`` API), buffer them in memory, and
emit one Parquet row group per batch. Centralising the implementation here
keeps the two callers in lockstep — fix a bug once, fix it everywhere.
"""

from __future__ import annotations

import typing as ty
from contextlib import contextmanager
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

# Default rows-per-row-group. Trades memory for row-group overhead; matches
# the value originally chosen for ``ParquetEntryWriter``.
DEFAULT_BATCH_SIZE = 64_000


class ParquetTable:
    """
    Minimal csv.writer-compatible facade around a streaming Parquet writer.

    Rows are buffered as tuples; every ``batch_size`` rows we transpose to
    per-column lists, build a ``pa.RecordBatch`` against the declared schema,
    and append it as a Parquet row group. Memory is O(batch_size * row_width).

    Exposes ``writerow`` / ``writerows`` so it drops into code paths that
    previously expected a ``csv.writer``.
    """

    def __init__(
        self, writer: pq.ParquetWriter, schema: pa.Schema, batch_size: int
    ) -> None:
        self._writer = writer
        self._schema = schema
        self._batch_size = batch_size
        self._buffer: ty.List[ty.Sequence] = []

    def writerow(self, row: ty.Sequence) -> None:
        self._buffer.append(tuple(row))
        if len(self._buffer) >= self._batch_size:
            self._flush()

    def writerows(self, rows: ty.Iterable[ty.Sequence]) -> None:
        buf = self._buffer
        batch_size = self._batch_size
        for row in rows:
            buf.append(tuple(row))
            if len(buf) >= batch_size:
                self._flush()

    def _flush(self) -> None:
        if not self._buffer:
            return
        columns = list(zip(*self._buffer))
        arrays = [
            pa.array(list(col), type=field.type)
            for col, field in zip(columns, self._schema)
        ]
        batch = pa.RecordBatch.from_arrays(arrays, schema=self._schema)
        self._writer.write_batch(batch)
        self._buffer.clear()

    def close(self) -> None:
        self._flush()
        self._writer.close()


@contextmanager
def parquet_writer(
    path: Path,
    schema: pa.Schema,
    batch_size: int = DEFAULT_BATCH_SIZE,
    compression: str = "zstd",
) -> ty.Iterator[ParquetTable]:
    """
    Open a single streaming Parquet file at ``path`` and yield a
    :class:`ParquetTable` bound to it. The file is flushed and closed on exit
    even if an exception propagates.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq_writer = pq.ParquetWriter(path, schema, compression=compression)
    table = ParquetTable(pq_writer, schema, batch_size)
    try:
        yield table
    finally:
        table.close()
