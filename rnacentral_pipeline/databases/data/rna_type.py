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
from attr.validators import optional
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.sequence_ontology import tree


@attr.s(frozen=True, slots=True)
class SoTermInfo:
    so_id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def ncRNA(cls):
        return cls("ncRNA", "SO:0000655")

    def is_a(self, value: str):
        return self.name == value or self.so_id == value


@attr.s(frozen=True, slots=True, hash=True)
class RnaType:
    so_term = attr.ib(validator=is_a(SoTermInfo))
    insdc = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def from_so_term(cls, so_tree, so_id) -> "RnaType":
        so_name = so_tree.name_to_id.get(so_id)
        so_node = so_tree[so_id]
        insdc = next(tree.insdc_synonyms(so_node), None)
        return cls(so_term=SoTermInfo(so_id=so_id, name=so_name), insdc=insdc)

    @classmethod
    def from_so_id(cls, so_tree, so_id) -> "RnaType":
        """
        Build and RNA type given the SO id, like SO:0000001
        """
        if so_id.startswith('SO:'):
            return cls.from_so_term(so_tree, so_id)
        return cls.from_so_term(so_tree, so_tree.name_to_id[so_id])

    @classmethod
    def from_insdc_term(cls, so_tree, term) -> "RnaType":
        """
        Build an RNA type given an INSDC RNA type.
        """
        if term == 'other' or term == 'ncRNA' or term == 'sRNA':
            return cls.from_so_term(so_tree, "SO:0000655")
        if term == 'misc_RNA':
            return cls.from_so_term(so_tree, "SO:0000673")
        name = None
        if term in so_tree.name_to_id:
            name = so_tree.name_to_id[term]
        elif term in so_tree.insdc_to_id:
            name = so_tree.insdc_to_id[term]
        else:
            raise ValueError(f"Unknown INSDC term {term}")
        return cls.from_so_term(so_tree, name)

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
