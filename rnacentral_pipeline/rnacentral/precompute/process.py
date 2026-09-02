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

import typing as ty
from contextlib import contextmanager
from pathlib import Path

import attr

from rnacentral_pipeline import psql, schemas, writers
from rnacentral_pipeline.parquet_writers import TypedParquetWrapper
from rnacentral_pipeline.rnacentral.precompute.data.context import Context
from rnacentral_pipeline.rnacentral.precompute.data.sequence import Sequence
from rnacentral_pipeline.rnacentral.precompute.data.update import SequenceUpdate

AnUpdate = SequenceUpdate


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
# str -> typed conversion before the row hits ``ParquetTable.writerow``. That
# bridging is the shared :class:`TypedParquetWrapper`; we wrap the underlying
# parquet table rather than touching the writeable producers (the CSV path
# stays bit-identical that way).

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
            precompute=TypedParquetWrapper(raw.precompute, schemas.PRECOMPUTE_DATA),
            qa=TypedParquetWrapper(raw.qa, schemas.PRECOMPUTE_QA),
        )


def parse(context_path: Path, data_path: Path) -> ty.Iterable[AnUpdate]:
    """
    Parse the given json file (handle) using the repeat tree at `repeat_path`,
    and produce an iterable of updates for the database.
    """

    context = Context.from_directory(context_path)
    with data_path.open("r") as handle:
        raw = psql.json_handler(handle)
        for sequence in raw:
            sequence = Sequence.build(context.so_tree, sequence)
            yield SequenceUpdate.from_sequence(context, sequence)
