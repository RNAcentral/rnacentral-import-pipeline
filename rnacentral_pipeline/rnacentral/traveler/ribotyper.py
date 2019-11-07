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

import os
import re

from pathlib import Path

import attr
from attr.validators import instance_of as is_a


class UnknownStrandName(Exception):
    pass


def as_strand(raw):
    if raw == 'plus':
        return 1
    if raw == 'minus':
        return -1
    raise UnknownStrandName(raw)


@attr.s()
class Result(object):
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

    @classmethod
    def from_result(cls, row):
        parts = re.split(r'\s+', row, maxsplit=24)
        if parts[2] == 'FAIL':
            return None
        return cls(
            target=parts[1],
            status=parts[2],
            length=parts[3],
            fm=parts[4],
            fam=parts[5],
            domain=parts[6],
            model=parts[7],
            strand=as_strand(parts[8]),
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
        )


def parse(filename):
    with open(filename, 'r') as raw:
        for line in raw:
            if line.startswith('#'):
                continue
            result = Result.from_result(line)
            if result:
                yield result


def as_dict(directory: Path):
    ribotyper_fn = '.ribotyper.long.out'
    basename = directory.name
    filenames = [Path(basename + ribotyper_fn), Path(ribotyper_fn)]
    for fn in filenames:
        filename = directory / fn
        if filename.exists():
            break
    else:
        raise ValueError("No fribotyper result file in: %s " % directory)
    return {p.target: p for p in parse(filename)}
