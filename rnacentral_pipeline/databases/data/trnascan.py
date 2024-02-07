# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from __future__ import annotations

import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


def intron_endpoint(raw: str) -> ty.Optional[int]:
    """Computes the endpoint of an intron. tRNAScan-SE uses '0' to mean no end
    point.

    >>> intron_endpoint("1")
    1
    >>> intron_endpoint("0")
    None
    """
    if raw == "0":
        return None
    return int(raw)


@attr.s()
class TRnaScanResults:
    """This represents a single result from tRNAScan-SE. This is essentially a
    literal translation of the results into a python object, with few extra
    features. The only thing to note is that the pseudogenes
    """

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
    def from_line(cls, line: str) -> TRnaScanResults:
        """Parse a given line from tRNAScan-SE results into a result object.

        >>> TRnaScanResults.from_line("URS0000CDDC17 	1	1 	91	Ser	TGA	0	0	26.8	pseudo")
        TRnaScanResults(sequence_id="URS0000C7FBE7", hit_index=1, sequence_start=1, sequence_stop=73, trna_type="Val", anticodon="TAC", intron_start=None, intron_stop=None, score=40.6, note="pseudo",)
        """
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
        """
        Check if this is marked as a pseudogene.
        """
        return "pseudo" in self.note
