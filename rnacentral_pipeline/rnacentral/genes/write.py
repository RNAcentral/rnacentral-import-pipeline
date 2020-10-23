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


def write_csv(rows, output, data_type: data.DataType, **kwargs) -> bool:
    out = output / f"{data_type.name}.csv"
    with out.open("a") as out:
        writer = csv.writer(out)
        rows = list(it.chain.from_iterable(rows))
        writer.writerows(rows)
    return True


def write_bed(bed, output: Path, data_type: data.DataType, extended=False, **kwargs) -> bool:
    out = output / f"{data_type.name}.bed"
    with out.open("a") as out:
        write_bed_text(it.chain.from_iterable(bed), out, extended=extended)
    return True


def write_gff(gff, output: Path, data_type: data.DataType, header=True) -> bool:
    out = output / f"{data_type.name}.gff"
    with out.open("a") as out:
        return write_gff_text(it.chain.from_iterable(gff), out, header=header, allow_no_features=True)


def write(
    results: ty.Iterable[data.FinalizedState],
    format: Format,
    path: Path,
    include_genes=True,
    allowed_members=None,
    allowed_data_types=data.DataType.all(),
    extended_bed=False,
):
    first = True
    written = False
    for result in results:
        for (name, locations) in result.data_types():
            if name not in allowed_data_types:
                LOGGER.debug("Skipping %s/%s since it is ignored", result.chromosome, name.name)
                continue

            if not data:
                LOGGER.debug("No entries in %s/%s", result.chromosome, name.name)
                continue

            method = None
            writer = None
            if format == Format.Csv:
                kwargs = {}
                if name == data.DataType.clustered:
                    kwargs['allowed_members'] = allowed_members
                method = op.methodcaller('as_writeable', **kwargs)
                writer = write_csv

            elif format == Format.Bed:
                kwargs = {}
                if name == data.DataType.clustered:
                    kwargs = {'allowed_members': allowed_members, 'include_gene': include_genes}
                method = op.methodcaller("as_bed", **kwargs)
                writer = partial(write_bed, extended=extended_bed)

            elif format == Format.Gff:
                kwargs = {}
                if name == data.DataType.clustered:
                    kwargs = {'allowed_members': allowed_members, 'include_gene': include_genes}
                method = op.methodcaller("as_features", **kwargs)
                writer = partial(write_gff, header=first)
            else:
                raise ValueError("Cannot write to format: %s" % format)

            entries = map(method, locations)
            data_written = writer(entries, path, name)
            written = written or data_written
            first = False

    if not written:
        raise ValueError("No features of any type written")
    return True
