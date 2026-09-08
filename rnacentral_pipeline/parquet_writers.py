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

import csv
import typing as ty
from contextlib import contextmanager
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from rnacentral_pipeline.output_format import is_parquet

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
        self._rows = 0

    def writerow(self, row: ty.Sequence) -> None:
        if len(row) != len(self._schema):
            raise ValueError(
                f"Row length ({len(row)}) does not match schema length ({len(self._schema)})"
            )
        self._buffer.append(tuple(row))
        self._rows += 1
        if len(self._buffer) >= self._batch_size:
            self._flush()

    def writerows(self, rows: ty.Iterable[ty.Sequence]) -> None:
        buf = self._buffer
        batch_size = self._batch_size
        schema_len = len(self._schema)
        for row in rows:
            if len(row) != schema_len:
                raise ValueError(
                    f"Row length ({len(row)}) does not match schema length ({len(self._schema)})"
                )
            buf.append(tuple(row))
            self._rows += 1
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

    @property
    def rows(self) -> int:
        """Rows written so far, buffered or flushed."""
        return self._rows

    def close(self) -> None:
        self._flush()
        self._writer.close()


def _to_int_or_none(value):
    if value in ("", "NaN", "nan", None):
        return None
    return int(value)


def _to_float_or_none(value):
    if value in ("", "NaN", "nan", None):
        return None
    return float(value)


def _to_bool(value):
    if value in (True, "1", "t", "T", "true", "True"):
        return True
    if value in (False, "0", "f", "F", "false", "False"):
        return False
    raise ValueError(f"Cannot interpret {value!r} as bool")


def _to_bool_or_none(value):
    if value in ("", None):
        return None
    return _to_bool(value)


def _to_string_list(value):
    """
    Accept either a real list or the Postgres array literal the CSV-era
    ``writeable()`` methods emit ('{}', '{"a","b"}') and return a list.

    pgloader parsed those literals into arrays on its way into text[] columns.
    DuckDB does not: it refuses to cast VARCHAR '{}' to VARCHAR[], so the
    parquet file has to carry a real list. Parsing here keeps that knowledge
    in one place instead of branching on the output format inside every
    writeable() that has an array column.
    """
    if value is None or value == "":
        return None
    if isinstance(value, (list, tuple)):
        return list(value)
    if not (value.startswith("{") and value.endswith("}")):
        raise ValueError(f"Cannot interpret {value!r} as an array")
    inner = value[1:-1]
    if not inner:
        return []
    # csv handles the quoting rules of the literal ('{"a,b","c"}') for us.
    return next(csv.reader([inner]))


def converter_for(field: pa.Field) -> ty.Callable[[ty.Any], ty.Any]:
    """
    Return a callable that converts a CSV-string row value into the typed
    value required by ``field.type``. Used by adapters that wrap a
    :class:`ParquetTable` to bridge the legacy stringly-typed ``writeable()``
    methods to typed Arrow batches.
    """
    t = field.type
    if pa.types.is_boolean(t):
        return _to_bool_or_none if field.nullable else _to_bool
    if pa.types.is_integer(t):
        return _to_int_or_none if field.nullable else (lambda v: int(v))
    if pa.types.is_floating(t):
        return _to_float_or_none if field.nullable else (lambda v: float(v))
    if pa.types.is_list(t):
        return (
            _to_string_list if field.nullable else (lambda v: _to_string_list(v) or [])
        )
    return lambda v: v


class TypedParquetWrapper:
    """
    Wrap a :class:`ParquetTable` with per-column converters so callers that
    emit string rows (from CSV-era ``writeable()`` methods) can write to a
    typed Parquet schema without changing the producer.
    """

    def __init__(self, table: "ParquetTable", schema: pa.Schema) -> None:
        self._table = table
        self._converters = [converter_for(f) for f in schema]

    def writerow(self, row: ty.Sequence) -> None:
        if len(row) != len(self._converters):
            raise ValueError(
                f"Row length ({len(row)}) does not match schema length ({len(self._converters)})"
            )
        self._table.writerow(tuple(c(v) for c, v in zip(self._converters, row)))

    def writerows(self, rows: ty.Iterable[ty.Sequence]) -> None:
        for row in rows:
            self.writerow(row)


@contextmanager
def typed_parquet_writer(
    path: Path,
    schema: pa.Schema,
    batch_size: int = DEFAULT_BATCH_SIZE,
    compression: str = "zstd",
) -> ty.Iterator[TypedParquetWrapper]:
    """
    Open a streaming parquet file and wrap it in :class:`TypedParquetWrapper`,
    so writers that emit stringly-typed rows can target a typed schema.
    """
    with parquet_writer(path, schema, batch_size, compression) as table:
        yield TypedParquetWrapper(table, schema)


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


@contextmanager
def row_writer(
    path: Path,
    schema: pa.Schema,
    csv_options: ty.Optional[ty.Mapping[str, ty.Any]] = None,
) -> ty.Iterator[ty.Any]:
    """
    Open ``path`` as either a CSV or a Parquet file, per the shared
    ``RNAC_OUTPUT_FORMAT`` switch, and yield a writer exposing the
    ``csv.writer`` API (``writerow``/``writerows``) in both cases.

    For the many small metadata parsers (ensembl assemblies, rfam clans, ols
    terms, ...) that emit stringly-typed rows straight from a query result,
    this is the whole format branch: ``with row_writer(path, SCHEMA) as w``.
    """
    if is_parquet():
        with typed_parquet_writer(path, schema) as writer:
            yield writer
        return

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        yield csv.writer(handle, **(csv_options or {}))
