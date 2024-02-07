# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


def intron_endpoint(raw: str) -> ty.Optional[int]:
    if raw == "0":
        return None
    return int(raw)


@attr.s()
class TRnaScanResults:
    sequence_id: str = attr.ib(validator=is_a(str))
    hit_index: int = attr.ib(validator=is_a(int))
    sequence_start: int = attr.ib(validator=is_a(int))
    sequence_stop: int = attr.ib(validator=is_a(int))
    anticodon: str = attr.ib(validator=is_a(str))
    trna_type: str = attr.ib(validator=is_a(str))
    intron_start: ty.Optional[int] = attr.ib(validator=optional(is_a(int)))
    intron_stop: ty.Optional[int] = attr.ib(validator=optional(is_a(int)))
    score: float = attr.ib(validator=is_a(float))
    note: str = attr.ib(validator=is_a(str))

    @classmethod
    def from_line(cls, line: str) -> "TRnaScanResults":
        parts = [p.strip() for p in line.split("\t")]
        return cls(
            sequence_id=parts[0],
            hit_index=int(parts[1]),
            sequence_start=int(parts[2]),
            sequence_stop=int(parts[3]),
            anticodon=parts[5],
            trna_type=parts[4],
            intron_start=intron_endpoint(parts[6]),
            intron_stop=intron_endpoint(parts[7]),
            score=float(parts[8]),
            note=parts[9],
        )

    @property
    def is_pseduo(self) -> bool:
        return "pseudo" in self.note
