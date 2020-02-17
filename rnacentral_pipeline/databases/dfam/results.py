# -*- coding: utf-8 -*-

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

import re
import csv
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.rfam.infernal_results import convert_strand


@attr.s()
class Result(object):
    model_name = attr.ib(validator=is_a(str))
    model_accession = attr.ib(validator=is_a(str))
    upi = attr.ib(validator=is_a(str))
    bits = attr.ib(validator=is_a(float))
    e_value = attr.ib(validator=is_a(float))
    bias = attr.ib(validator=is_a(float))
    model_start = attr.ib(validator=is_a(int))
    model_end = attr.ib(validator=is_a(int))
    strand = attr.ib(validator=is_a(int))
    sequence_start = attr.ib(validator=is_a(int))
    sequence_end = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, parts):
        return cls(
            model_name=parts[0],
            model_accession=parts[1],
            upi=parts[2],
            bits=float(parts[3]),
            e_value=float(parts[4]),
            bias=float(parts[5]),
            model_start=int(parts[6]),
            model_end=int(parts[7]),
            strand=convert_strand(parts[8]),
            sequence_start=int(parts[9]),
            sequence_end=int(parts[10]),
        )

    def writeable(self):
        return [
            self.upi,
            self.model_accession,
            self.sequence_start,
            self.sequence_end,
            self.strand,
            self.model_start,
            self.model_end,
            self.e_value,
            self.bits,
        ]


def parse(handle):
    for line in handle:
        if line.startswith('#'):
            continue
        parts = re.split(r'\s+', line.strip(), 15)
        yield Result.build(parts)


def write(handle, output):
    writer = csv.writer(output)
    data = parse(handle)
    data = map(op.methodcaller('writeable'), data)
    writer.writerows(data)
