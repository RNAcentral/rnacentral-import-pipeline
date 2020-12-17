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
from functools import lru_cache

import attr
import networkx as nx

from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree

NORMALIZED_IDS = {
    "SO:0000253",  # tRNA
    "SO:0000274",  # snRNA
    "SO:0000275",  # snoRNA
    "SO:0000375",  # 5.8S rRNA
    "SO:0000650",  # SSU rRNA
    "SO:0000651",  # LSU rRNA
    "SO:0000652",  # 5S rRNA
    "SO:0002128",  # mito rRNA
}


@attr.s(auto_attribs=True, frozen=True, slots=True)
class SoTermInfo:
    name: str
    so_id: str

    @classmethod
    def ncRNA(cls):
        return cls("ncRNA", "SO:0000655")

    def is_a(self, value: str):
        return self.name == value or self.so_id == value


@attr.s(auto_attribs=True, frozen=True, slots=True)
class NormalizedSoTermInfo:
    so_term: SoTermInfo
    is_normalized: bool

    def __getattr__(self, name):
        return getattr(self.so_term, name)


SoTree = ty.Tuple[SoTermInfo]


@lru_cache()
def normalized_term(so_id: str, ontology) -> NormalizedSoTermInfo:
    for parent_so_id in NORMALIZED_IDS:
        if parent_so_id == so_id:
            so_term = SoTermInfo(so_id, ontology.id_to_name[so_id])
            return NormalizedSoTermInfo(so_term, True)
        paths = nx.all_simple_paths(ontology, source=so_id, target=parent_so_id)
        paths = list(paths)
        if len(paths) == 1:
            name = ontology.id_to_name[parent_so_id]
            so_term = SoTermInfo(name, parent_so_id)
            return NormalizedSoTermInfo(so_term, True)

    name = ontology.id_to_name[so_id]
    so_term = SoTermInfo(name, so_id)
    return NormalizedSoTermInfo(so_term, False)


@attr.s(auto_attribs=True, frozen=True, slots=True, hash=True)
class RnaType:
    insdc: str
    so_term: str
    ontology_terms: SoTree
    normalized_term: SoTermInfo

    @classmethod
    def build(cls, insdc, so_term, ontology):
        subtree = so_tree.rna_type_tree(ontology, so_term)
        subtree = tuple([SoTermInfo(p[1], p[0]) for p in subtree])
        so_id = ontology.name_to_id[so_term]
        return cls(
            insdc=insdc,
            so_term=so_term,
            ontology_terms=subtree,
            normalized_term=normalized_term(so_id, ontology),
        )

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


@attr.s(auto_attribs=True, frozen=True, slots=True, hash=True)
class RnaTypeGraph:
    graph: nx.Graph

    @classmethod
    def from_containers(cls, containers):
        graph = nx.Graph()
        for container in containers:
            rna_type = container.rna_type
            node_name = rna_type.so_term.name
            graph.add_node(node_name, {"data", []})
            graph[node_name]["data"].append(container)
            for so_term in rna_type.ontology_terms:
                graph.add_node(so_term.name)
            ns = [n.name for n in rna_type.ontology_terms]
            if len(ns) > 1:
                edges = tuple(zip(ns, ns[1:]))
                graph.add_edges_from(edges)

        cls(graph=graph)
