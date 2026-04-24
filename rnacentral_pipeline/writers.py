# -*- coding: utf-8 -*-

from __future__ import annotations

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import csv
import logging
import os
import typing as ty
from contextlib import ExitStack, contextmanager
from pathlib import Path

import attr
import pyarrow as pa
import pyarrow.parquet as pq

from rnacentral_pipeline import schemas
from rnacentral_pipeline.databases import data

_C = ty.TypeVar("_C")

LOGGER = logging.getLogger(__name__)

# Number of rows buffered per logical table before we emit a Parquet row group.
# Chosen as a middle ground between memory footprint and row-group overhead;
# profile during the spike and retune if needed.
PARQUET_ROW_GROUP_SIZE = 64_000


@attr.s()
class EntryWriter:
    accessions = attr.ib()
    short_sequences = attr.ib(
        metadata={"csv_options": {"delimiter": ",", "lineterminator": "\n"}}
    )
    long_sequences = attr.ib(
        metadata={"csv_options": {"delimiter": ",", "lineterminator": "\n"}}
    )
    references = attr.ib()
    ref_ids = attr.ib()
    regions = attr.ib()
    secondary_structure = attr.ib()
    related_sequences = attr.ib()
    features = attr.ib()
    interactions = attr.ib()
    terms = attr.ib()
    go_annotations = attr.ib()
    go_publication_mappings = attr.ib()

    def write(self, entries: ty.Iterable[data.Entry]):
        total = 0
        invalid = 0
        for entry in entries:
            total += 1
            if not entry.is_valid():
                invalid += 1
                continue

            self.accessions.writerows(entry.write_ac_info())
            self.short_sequences.writerows(entry.write_seq_short())
            self.long_sequences.writerows(entry.write_seq_long())
            self.references.writerows(entry.write_refs())
            self.ref_ids.writerows(entry.write_ref_ids())
            self.regions.writerows(entry.write_sequence_regions())
            self.secondary_structure.writerows(entry.write_secondary_structure())
            self.related_sequences.writerows(entry.write_related_sequences())
            self.features.writerows(entry.write_sequence_features())
            self.interactions.writerows(entry.write_interactions())
            self.terms.writerows(entry.write_ontology_terms())

            for annotations in entry.go_annotations:
                self.go_annotations.writerows(annotations.writeable())
                self.go_publication_mappings.writerows(
                    annotations.writeable_publication_mappings()
                )
                self.terms.writerows(annotations.writeable_ontology_terms())

        if not total:
            raise ValueError("Found no entries to write")

        LOGGER.info("Wrote %i entries", total)
        if invalid:
            LOGGER.warn("Did not write %i of %i total entries", invalid, total)
            if invalid == total:
                raise ValueError("No valid entries to write")


@attr.s()
class OntologyAnnnotationWriter:
    go_annotations = attr.ib()
    go_publication_mappings = attr.ib()
    terms = attr.ib()

    def write(self, annotations: ty.Iterable[data.GoTermAnnotation]):
        for anno in annotations:
            self.go_annotations.writerows(anno.writeable())
            self.go_publication_mappings.writerows(
                anno.writeable_publication_mappings()
            )
            self.terms.writerows(anno.writeable_ontology_terms())


@contextmanager
def build(cls: ty.Type[_C], path: Path) -> ty.Iterator[_C]:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    handles = {}
    with ExitStack() as stack:
        for field in attr.fields(cls):
            out = path / f"{field.name}.csv"
            handle = Path(out).open("w")
            options = field.metadata.get(
                "csv_options",
                {
                    "delimiter": ",",
                    "quotechar": '"',
                    "quoting": csv.QUOTE_ALL,
                    "lineterminator": "\n",
                },
            )
            handles[field.name] = csv.writer(handle, **options)
            stack.enter_context(handle)
        yield cls(**handles)


def entry_writer(path: Path):
    """
    Open the standard entry writer at ``path``. Honours the
    ``RNAC_OUTPUT_FORMAT`` environment variable: set to ``parquet`` to emit
    streaming Parquet files (see :func:`parquet_entry_writer`). Defaults to
    the historical CSV path.
    """
    if os.environ.get("RNAC_OUTPUT_FORMAT", "csv").lower() == "parquet":
        return parquet_entry_writer(path)
    return build(EntryWriter, path)


# ---------------------------------------------------------------------------
# Parquet path
# ---------------------------------------------------------------------------


class _ParquetTable:
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


@attr.s()
class ParquetEntryWriter:
    """
    Drop-in API-compatible replacement for :class:`EntryWriter` that emits one
    streaming Parquet file per logical table instead of 13 CSVs. Callers use
    :func:`parquet_entry_writer` exactly like :func:`entry_writer`.
    """

    accessions = attr.ib()
    short_sequences = attr.ib()
    long_sequences = attr.ib()
    references = attr.ib()
    ref_ids = attr.ib()
    regions = attr.ib()
    secondary_structure = attr.ib()
    related_sequences = attr.ib()
    features = attr.ib()
    interactions = attr.ib()
    terms = attr.ib()
    go_annotations = attr.ib()
    go_publication_mappings = attr.ib()

    def write(self, entries: ty.Iterable[data.Entry]):
        total = 0
        invalid = 0
        for entry in entries:
            total += 1
            if not entry.is_valid():
                invalid += 1
                continue

            self.accessions.writerows(entry.write_ac_info())
            self.short_sequences.writerows(entry.write_seq_short())
            self.long_sequences.writerows(entry.write_seq_long())
            self.references.writerows(entry.write_refs())
            self.ref_ids.writerows(entry.write_ref_ids())
            self.regions.writerows(entry.write_sequence_regions())
            self.secondary_structure.writerows(entry.write_secondary_structure())
            self.related_sequences.writerows(entry.write_related_sequences())
            self.features.writerows(entry.write_sequence_features())
            self.interactions.writerows(entry.write_interactions())
            self.terms.writerows(entry.write_ontology_terms())

            for annotations in entry.go_annotations:
                self.go_annotations.writerows(annotations.writeable())
                self.go_publication_mappings.writerows(
                    annotations.writeable_publication_mappings()
                )
                self.terms.writerows(annotations.writeable_ontology_terms())

        if not total:
            raise ValueError("Found no entries to write")

        LOGGER.info("Wrote %i entries", total)
        if invalid:
            LOGGER.warn("Did not write %i of %i total entries", invalid, total)
            if invalid == total:
                raise ValueError("No valid entries to write")


@contextmanager
def parquet_entry_writer(
    path: Path,
    batch_size: int = PARQUET_ROW_GROUP_SIZE,
    compression: str = "zstd",
) -> ty.Iterator[ParquetEntryWriter]:
    """
    Open one streaming Parquet writer per logical ``EntryWriter`` table and
    yield a ``ParquetEntryWriter`` bound to them. On exit, each writer is
    flushed and closed even if an exception propagates.
    """

    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    tables: ty.Dict[str, _ParquetTable] = {}
    with ExitStack() as stack:
        for name, schema in schemas.ENTRY_WRITER_SCHEMAS.items():
            out = path / f"{name}.parquet"
            pq_writer = pq.ParquetWriter(out, schema, compression=compression)
            table = _ParquetTable(pq_writer, schema, batch_size)
            stack.callback(table.close)
            tables[name] = table
        yield ParquetEntryWriter(**tables)
