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

import enum
import itertools as it

from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import \
    write_bed_text
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.gff3 import \
    write_gff_text


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


def write_csv(data, output):
    rows = it.chain.from_iterable(d.writeable() for d in data)
    writer = csv.writer(output)
    writer.writerows(rows)


def write_bed(data, output, extended=False):
    bed = it.chain.from_iterable(d.as_bed() for d in data)
    write_bed_text(bed, output, extended=extended)


def write_gff(data, output):
    gff = it.chain.from_iterable(d.as_features() for d in data)
    write_gff_text(gff, output)


def write(data, format, output):
    if format == Format.Csv:
        return write_csv(data, output)
    elif format == Format.Bed:
        return write_bed(data, output)
    elif format == Format.Gff:
        return write_gff(data, output)
    raise ValueError("Cannot write to format: %s" % format)
