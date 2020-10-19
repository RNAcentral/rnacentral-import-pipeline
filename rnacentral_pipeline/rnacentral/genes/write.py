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


def write_csv(rows, output: Path, name: str):
    out = output / f"{name}.csv"
    with out.open("w") as out:
        writer = csv.writer(out)
        writer.writerows(rows)


def write_bed(bed, output: Path, name: str, extended=False):
    out = output / f"{name}.bed"
    with out.open("w") as out:
        write_bed_text(bed, out, extended=extended)


def write_gff(gff, output: Path, name: str):
    out = output / f"{name}.gff"
    with out.open("w") as out:
        return write_gff_text(gff, out, allow_no_features=True)


def apply_handlers(
    finished: ty.Iterable[data.Finalized], locus_handler, rejected_handler, path, writer
):
    written = False
    for finalized in finished:
        for locus in finalized.locuses:
            result = locus_handler(locus)
            locus_written = writer(result, path, "genes")
            writen = written or locus_written

        for rejected in finalized.rejected:
            result = rejected_handler(rejected)
            rejected_written = writer(result, path, "rejected")
            written = written or rejected_written
    return written


def write(
    data: ty.Iterable[data.Finalized],
    format: Format,
    path: Path,
    include_genes=None,
    include_representative=None,
    include_members=None,
    extended_bed=None,
    include_rejected=False,
):
    if format == Format.Csv:
        locus_method = op.methodcaller(
            "as_writeable",
            include_representative=include_representative,
            include_members=include_members,
        )
        rejected_method = op.methodcaller("as_writeable")
        apply_handlers(data, locus_method, rejected_method, path, write_csv)
        return True

    elif format == Format.Bed:
        locus_method = op.methodcaller(
            "as_bed",
            include_gene=include_genes,
            include_representative=include_representative,
            include_members=include_members,
        )
        rejected_method = op.methodcaller("as_bed")
        writer = partial(write_bed, extended=extended_bed)
        apply_handlers(data, locus_method, rejected_method, path, writer)
        return True

    elif format == Format.Gff:
        locus_method = op.methodcaller(
            "as_features",
            include_gene=include_genes,
            include_representative=include_representative,
            include_members=include_members,
        )
        rejected_method = op.methodcaller("as_features")
        written = apply_handlers(data, locus_method, rejected_method, path, write_gff)
        if not written:
            raise ValueError("No features of any time written")
        return True

    raise ValueError("Cannot write to format: %s" % format)
