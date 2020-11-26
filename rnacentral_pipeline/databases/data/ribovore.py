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
from attr.validators import instance_of as is_a

from .regions import UnknownStrand


@attr.s()
class RibovoreResult(object):
    target: str = attr.ib(validator=is_a(str))
    status: str = attr.ib(validator=is_a(str))
    length: int = attr.ib(validator=is_a(int), converter=int)
    fm: int = attr.ib(validator=is_a(int), converter=int)
    fam: str = attr.ib(validator=is_a(str))
    domain: str = attr.ib(validator=is_a(str))
    model: str = attr.ib(validator=is_a(str))
    strand: int = attr.ib(validator=is_a(int))
    ht: int = attr.ib(validator=is_a(int), converter=int)
    tscore: float = attr.ib(validator=is_a(float), converter=float)
    bscore: float = attr.ib(validator=is_a(float), converter=float)
    bevalue: float = attr.ib(validator=is_a(float), converter=float)
    tcov: float = attr.ib(validator=is_a(float), converter=float)
    bcov: float = attr.ib(validator=is_a(float), converter=float)
    bfrom: int = attr.ib(validator=is_a(int), converter=int)
    bto: int = attr.ib(validator=is_a(int), converter=int)
    mfrom: int = attr.ib(validator=is_a(int), converter=int)
    mto: int = attr.ib(validator=is_a(int), converter=int)
    model_length = attr.ib(validator=optional(is_a(int)), default=None)

    @classmethod
    def from_result_line(cls, row: str, lengths=None) -> ty.Optional["RibovoreResult"]:
        parts = re.split(r'\s+', row, maxsplit=24)
        if parts[2] == 'FAIL':
            return None
        strand = None
        if parts[8] == 'plus':
            strand = 1
        elif parts[8] == 'minus':
            strand = -1
        else:
            raise UnknownStrand(parts[8])

        model_length = None
        if lengths:
            model_length = lengths.get(parts[7], None)

        return cls(
            target=parts[1],
            status=parts[2],
            length=parts[3],
            fm=parts[4],
            fam=parts[5],
            domain=parts[6],
            model=parts[7],
            strand=strand,
            ht=parts[9],
            tscore=parts[10],
            bscore=parts[11],
            bevalue=parts[13],
            tcov=parts[14],
            bcov=parts[15],
            bfrom=parts[16],
            bto=parts[17],
            mfrom=parts[18],
            mto=parts[19],
            model_length=model_length,
        )

    @property
    def model_coverage(self) -> ty.Optional[float]:
        if self.model_length:
            return float(self.mto - self.mfrom) / float(self.model_length)
        return None

    @property
    def sequence_coverage(self):
        return self.tcov
