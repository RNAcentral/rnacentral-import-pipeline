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

import itertools as it
import typing as ty

import attr
from attr.validators import instance_of as is_a


@attr.s(frozen=True, slots=True)
class SoTermInfo:
    name = attr.ib(validator=is_a(str))
    so_id = attr.ib(validator=is_a(str))

    @classmethod
    def ncRNA(cls):
        return cls("ncRNA", "SO:0000655")

    def is_a(self, value: str):
        return self.name == value or self.so_id == value


SoTree = ty.Tuple[SoTermInfo]


@attr.s(frozen=True, slots=True, hash=True)
class RnaType:
    insdc = attr.ib(validator=is_a(str))
    so_term = attr.ib(validator=is_a(str))
    ontology_terms = attr.ib(validator=is_a(tuple), type=SoTree)

    def is_a(self, so_name):
        return (
            self.insdc == so_name
            or self.so_term == so_name
            or any(p.is_a(so_name) for p in self.ontology_terms)
        )

    def common(self, other: "RnaType") -> ty.Optional[SoTermInfo]:
        zipped = zip(self.ontology_terms, other.ontology_terms)
        common = list(it.takewhile(lambda z: z[0] == z[1], zipped))
        if not common:
            return None
        return common

    def is_parent_of(self, other: "RnaType") -> bool:
        if len(other.ontology_terms) >= len(self.ontology_terms):
            return False

        pairs = zip(self.ontology_terms, other.ontology_terms)
        for (left, right) in pairs:
            if left != right:
                return False
        return True
