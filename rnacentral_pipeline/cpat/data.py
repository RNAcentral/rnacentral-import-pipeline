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
class CpatResult:
    urs_taxid = attr.ib(validator=is_a(str))
    fickett_score = attr.ib(validator=is_a(float))
    hexamer_score = attr.ib(validator=is_a(float))
    coding_prob = attr.ib(validator=is_a(float))
    protein_coding = attr.ib(validator=is_a(bool))

    @classmethod
    def build(cls, raw: ty.Dict[str, str], model_name: str, cutoffs: CpatCutoffs) -> CpatResult:
        prob = float(raw['Coding_prob'])
        coding = cutoffs.is_protein_coding(model_name, prob)
        return cls(
            urs_taxid=raw['seq_ID'],
            fickett_score=float(raw['Fickett']),
            hexamer_score=float(raw['Hexamer']),
            coding_prob=prob,
            protein_coding=coding,
        )

    def writeable(self) -> ty.List[str]:
        return [
            self.urs_taxid,
            str(self.fickett_score),
            str(self.hexamer_score),
            str(self.coding_prob),
            str(self.protein_coding),
        ]


@attr.s()
class CpatCutoffs:
    cutoffs: ty.Dict[str, float] = attr.ib(validator=is_a(dict), factory=dict)

    def add_cutoff(self, source: str, cutoff: float):
        self.cutoffs[source] = cutoff

    def is_protein_coding(self, source: str, coding_prob: float) -> bool:
        return coding_prob >= self.cutoffs[source]
