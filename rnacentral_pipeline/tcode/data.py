from __future__ import annotations

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


@attr.s()
class TcodeResult:
    urs = attr.ib(validator=is_a(str))
    length = attr.ib(validator=is_a(str))
    mean_score = attr.ib(validator=is_a(str))
    std_score = attr.ib(validator=is_a(str))
    is_protein_coding = attr.ib(validator=is_a(bool))

    @classmethod
    def build(
        cls,
        urs: str,
        size: ty.Optional[int],
        mean_score: ty.Optional[float],
        std_score: ty.Optional[float],
    ) -> TcodeResult:
        is_protein_coding = mean_score is not None and mean_score > 0.95
        return cls(
            urs=urs,
            length=_format_value(size),
            mean_score=_format_value(mean_score),
            std_score=_format_value(std_score),
            is_protein_coding=is_protein_coding,
        )

    def writeable(self) -> ty.List[str]:
        return [
            self.urs,
            self.length,
            self.mean_score,
            self.std_score,
            str(self.is_protein_coding),
        ]


@attr.s()
class TcodeWriter:
    results = attr.ib()

    def write(self, results: ty.Iterable[TcodeResult]):
        for result in results:
            self.results.writerow(result.writeable())


def _format_value(value: ty.Optional[ty.Union[int, float]]) -> str:
    if value is None:
        return ""
    return str(value)
