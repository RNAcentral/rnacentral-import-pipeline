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

import re
import typing as ty

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data.regions import UnknownStrand


def maybe(convert):
    def fn(value):
        if value is None or value == "-":
            return None
        return convert(value)

    return fn


maybe_int = maybe(int)
maybe_float = maybe(float)


@attr.s()
class RibovoreResult(object):
    target = attr.ib(validator=is_a(str))
    status = attr.ib(validator=is_a(str))
    length = attr.ib(validator=is_a(int))
    fm = attr.ib(validator=is_a(int))
    fam = attr.ib(validator=is_a(str))
    domain = attr.ib(validator=is_a(str))
    model = attr.ib(validator=is_a(str))
    strand = attr.ib(validator=optional(is_a(int)))
    ht = attr.ib(validator=optional(is_a(int)))
    tscore = attr.ib(validator=optional(is_a(float)))
    bscore = attr.ib(validator=optional(is_a(float)))
    bevalue = attr.ib(validator=optional(is_a(float)))
    tcov = attr.ib(validator=optional(is_a(float)))
    bcov = attr.ib(validator=optional(is_a(float)))
    bfrom = attr.ib(validator=optional(is_a(int)))
    bto = attr.ib(validator=optional(is_a(int)))
    mfrom = attr.ib(validator=optional(is_a(int)))
    mto = attr.ib(validator=optional(is_a(int)))
    model_length = attr.ib(validator=optional(is_a(int)), default=None)

    @classmethod
    def from_result_line(cls, row: str, lengths=None) -> "RibovoreResult":
        parts = re.split(r"\s+", row, maxsplit=24)
        strand = None
        if parts[8] == "plus":
            strand = 1
        elif parts[8] == "minus":
            strand = -1
        elif parts[8] == "-":
            strand = None
        else:
            raise UnknownStrand(parts[8])

        model_length = None
        if lengths:
            model_length = lengths.get(parts[7], None)

        return cls(
            target=parts[1],
            status=parts[2],
            length=int(parts[3]),
            fm=int(parts[4]),
            fam=parts[5],
            domain=parts[6],
            model=parts[7],
            strand=strand,
            ht=maybe_int(parts[9]),
            tscore=maybe_float(parts[10]),
            bscore=maybe_float(parts[11]),
            bevalue=maybe_float(parts[13]),
            tcov=maybe_float(parts[14]),
            bcov=maybe_float(parts[15]),
            bfrom=maybe_int(parts[16]),
            bto=maybe_int(parts[17]),
            mfrom=maybe_int(parts[18]),
            mto=maybe_int(parts[19]),
            model_length=model_length,
        )

    @property
    def model_coverage(self) -> ty.Optional[float]:
        length = self.modeled_length
        if length is not None:
            return float(length) / float(self.model_length)
        return None

    @property
    def modeled_length(self) -> ty.Optional[int]:
        if self.mto is not None and self.mfrom is not None
            return self.mto - self.mfrom
        return None

    @property
    def sequence_coverage(self):
        return self.tcov
