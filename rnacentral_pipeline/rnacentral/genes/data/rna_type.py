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

from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree


@attr.s(auto_attribs=True, frozen=True, slots=True)
class SoTermInfo:
    name: str
    so_id: str

    def is_a(self, value: str):
        return self.name == value or self.so_id == value


@attr.s(auto_attribs=True, frozen=True, slots=True)
class RnaType:
    insdc: str
    so_term: str
    ontology_terms: ty.List[SoTermInfo]

    def is_a(self, so_name):
        return (
            self.insdc == so_name
            or self.so_term == so_name
            or any(p.is_a(so_name) for p in self.ontology_terms)
        )

    def common(self, other: 'RnaType') -> ty.Optional[SoTermInfo]:
        zipped = zip(self.ontology_terms, other.ontology_terms)
        common = list(it.takewhile(lambda z: z[0] == z[1], zipped))
        if not common:
            return None
        return common
