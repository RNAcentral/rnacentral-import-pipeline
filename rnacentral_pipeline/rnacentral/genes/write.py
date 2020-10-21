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

import csv
import enum
import itertools as it
import operator as op
import typing as ty
from functools import partial
from pathlib import Path

from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import \
    write_bed_text
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.gff3 import \
    write_gff_text

from . import data


@enum.unique
class Format(enum.Enum):
    Csv = enum.auto()
    Bed = enum.auto()
    Gff = enum.auto()

    @classmethod
    def from_name(cls, name):
        for value in cls:
            if value.name.lower() == name:
                return value
        raise ValueError("Unknown Format %s" % name)

    @classmethod
    def names(cls):
        return [x.name for x in cls]


def write_csv(rows, output: Path, name: str) -> bool:
    out = output / f"{name}.csv"
    with out.open("w") as out:
        writer = csv.writer(out)
        first = next(rows, None)
        if first is None:
            return False
        writer.writerow(first)
        writer.writerows(rows)
    return True


def write_bed(bed, output: Path, name: str, extended=False) -> bool:
    out = output / f"{name}.bed"
    with out.open("w") as out:
        write_bed_text(bed, out, extended=extended)
    return True


def write_gff(gff, output: Path, name: str) -> bool:
    out = output / f"{name}.gff"
    with out.open("w") as out:
        return write_gff_text(gff, out, allow_no_features=True)


def apply_handlers(
    finished: ty.Iterable[data.FinalizedState],
    handlers,
    path: Path,
    writer,
    allowed={"genes", "rejected", "ignored"},
):
    written = False
    for finalized in finished:
        for (name, data) in finalized.data_types():
            if name not in allowed:
                continue
            if not data:
                continue
            handler = handlers[name]
            result = (handler(d) for d in data)
            data_written = writer(result, path, name)
            written = written or data_written
    return written


def write(
    data: ty.Iterable[data.FinalizedState],
    format: Format,
    path: Path,
    include_genes=None,
    allowed_members=None,
    allowed_data_types={"genes", "rejected", "ignored"},
    extended_bed=None,
):
    writer = None
    method_name = None
    if format == Format.Csv:
        method_name = "as_writeable"
        writer = write_csv

    elif format == Format.Bed:
        method_name = "as_bed"
        writer = partial(write_bed, extended=extended_bed)

    elif format == Format.Gff:
        method_name = "as_features"
        writer = write_gff

    else:
        raise ValueError("Cannot write to format: %s" % format)

    handlers = {
        "genes": op.methodcaller(
            method_name, include_gene=include_genes, allowed_members=allowed_members,
        ),
        "rejected": op.methodcaller(method_name),
        "ignored": op.methodcaller(method_name),
    }
    written = apply_handlers(data, handlers, path, writer, allowed=allowed_data_types)
    if not written:
        raise ValueError("No features of any type written")
    return True
