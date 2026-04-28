# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import itertools as it
import operator as op
import typing as ty
from contextlib import contextmanager
from pathlib import Path

import attr
import pyarrow as pa

from rnacentral_pipeline import psql, schemas, writers
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.data.update import (
    GenericUpdate,
    SequenceUpdate,
)

AnUpdate = ty.Union[SequenceUpdate, GenericUpdate]


@attr.s()
class Writer:
    precompute = attr.ib()
    qa = attr.ib()

    def write(self, updates: ty.Iterable[AnUpdate]):
        for update in updates:
            self.precompute.writerows(update.as_writeables())
            self.qa.writerows(update.writeable_statuses())


# ---------------------------------------------------------------------------
# Parquet path
#
# ``as_writeables`` / ``writeable_statuses`` emit string rows (ints and bools
# stringified for csv.writer). The parquet schemas in
# :mod:`rnacentral_pipeline.schemas` are typed, so each column needs a
# str -> typed conversion before the row hits ``ParquetTable.writerow``. We
# build one converter list per schema and wrap the underlying parquet table
# rather than touching the writeable producers (the CSV path stays
# bit-identical that way).


def _to_int_or_none(value: str) -> ty.Optional[int]:
    return int(value) if value not in ("", None) else None


def _to_bool_01(value: str) -> bool:
    # as_writeables stringifies ``int(bool_value)`` so we only ever see "0"
    # or "1" here. Defensive against pgloader-style "t"/"f" just in case.
    if value in ("1", "t", "true", "True"):
        return True
    if value in ("0", "f", "false", "False"):
        return False
    raise ValueError(f"Cannot interpret {value!r} as bool")


def _converter_for(field: pa.Field) -> ty.Callable[[str], ty.Any]:
    t = field.type
    if pa.types.is_boolean(t):
        return _to_bool_01
    if pa.types.is_integer(t):
        return _to_int_or_none if field.nullable else (lambda v: int(v))
    return lambda v: v


@attr.s()
class _TypedParquetWrapper:
    """
    Wraps a ``parquet_writers.ParquetTable`` with per-column converters so
    string rows from ``as_writeables`` / ``writeable_statuses`` line up with a
    typed parquet schema.
    """

    table = attr.ib()
    converters: ty.List[ty.Callable[[str], ty.Any]] = attr.ib()

    def writerow(self, row):
        self.table.writerow(tuple(c(v) for c, v in zip(self.converters, row)))

    def writerows(self, rows):
        for row in rows:
            self.writerow(row)


_FIELD_SCHEMAS = {
    "precompute": schemas.PRECOMPUTE_DATA,
    "qa": schemas.PRECOMPUTE_QA,
}


@contextmanager
def parquet_writer(path: Path) -> ty.Iterator[Writer]:
    """
    Open a precompute :class:`Writer` whose ``precompute`` and ``qa`` tables
    write to streaming Parquet files under ``path``.
    """
    with writers.build_parquet(Writer, path, _FIELD_SCHEMAS) as raw:
        yield Writer(
            precompute=_TypedParquetWrapper(
                raw.precompute,
                [_converter_for(f) for f in schemas.PRECOMPUTE_DATA],
            ),
            qa=_TypedParquetWrapper(
                raw.qa,
                [_converter_for(f) for f in schemas.PRECOMPUTE_QA],
            ),
        )


def parse(context_path: Path, data_path: Path) -> ty.Iterable[AnUpdate]:
    """
    Parse the given json file (handle) using the repeat tree at `repeat_path`,
    and produce an iterable of updates for the database.
    """

    context = Context.from_directory(context_path)
    with data_path.open("r") as handle:
        raw = psql.json_handler(handle)
        grouped = it.groupby(raw, op.itemgetter("upi"))
        for _, sequences in grouped:
            updates = []
            for sequence in sequences:
                sequence = Sequence.build(context.so_tree, sequence)
                update = SequenceUpdate.from_sequence(context, sequence)
                updates.append(update)
                yield update
            yield GenericUpdate.from_updates(context, updates)
