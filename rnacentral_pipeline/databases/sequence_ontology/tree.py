# -*- coding: utf-8 -*-

from __future__ import annotations

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
import logging
import typing as ty
from functools import lru_cache

import attr
from attr.validators import instance_of as is_a

from sqlitedict import SqliteDict
import networkx as nx
import obonet


REMOTE_ONTOLOGY = "https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so-simple.obo"

LOGGER = logging.getLogger(__name__)

BASE_SO_TERMS = [
    "ncRNA",
    "intron",
    "mRNA_region",
    "transcript",
]

RENAME = {
    ""
}

ALTERNATES = {
    'SO:0000380': [
        ("SO:0000673", "transcript"),
        ('SO:0000374', 'ribozyme'),
        ('SO:0000380', 'hammerhead_ribozyme'),
    ],
    'SO:0000209': [
        ('SO:0000655', 'ncRNA'),
        ('SO:0000252', 'rRNA'),
        ('SO:0000655', 'rRNA_primary_transcript'),
    ],
    'SO:0001904': [
        ("SO:0000655", "ncRNA"),
        ("SO:0001877", "lncRNA"),
        ("SO:0001904", "antisense_lncRNA"),
    ],
    'SO:0001244': [
        ("SO:0000673", "transcript"),
        ('SO:0001244', "pre_miRNA"),
    ]
}

SKIPPED_TERMS = {
    "SO:0000372",
    "SO:0000990",
    "SO:0000383",
    "SO:0000376",
    "SO:0000387",
    "SO:0000378",
    "SO:0000377",
    "SO:0000389",
    "SO:0000384",
    "SO:0000388",
    "SO:0000379",
}


@attr.s(frozen=True, hash=False)
class SoOntology:
    graph = attr.ib(validator=is_a(nx.MultiDiGraph))
    id_to_name: ty.Dict[str, ty.Optional[str]] = attr.ib(validator=is_a(dict))
    name_to_id: ty.Dict[str, str] = attr.ib(validator=is_a(dict))
    insdc_to_id: ty.Dict[str, str] = attr.ib(validator=is_a(dict))

    @classmethod
    def from_file(cls, filename) -> SoOntology:
        ont = obonet.read_obo(filename, ignore_obsolete=False)
        id_to_name = {}
        name_to_id = {}
        insdc_to_id = {}

        for old_id, old_node in ont.nodes(data=True):
            if "replaced_by" not in old_node:
                name = old_node.get("name", None)
                id_to_name[old_id] = name
                if name:
                    name_to_id[name] = old_id
                for insdc in insdc_synonyms(old_node):
                    insdc_to_id[insdc] = old_id
                continue

            for new_id in old_node["replaced_by"]:
                new_node = ont.nodes[new_id]
                if "name" in old_node:
                    name_to_id[old_node["name"]] = new_id

                id_to_name[old_id] = new_node.get("name")
                for insdc in insdc_synonyms(old_node):
                    insdc_to_id[insdc] = new_id

        return cls(graph=ont, id_to_name=id_to_name, name_to_id=name_to_id, insdc_to_id=insdc_to_id)

    def nodes(self):
        return self.graph.nodes(data=True)

    def as_node_id(self, id: str) -> str:
        if self.graph.has_node(id):
            if self.graph.nodes[id].get("replaced_by", []):
                return self.graph.nodes[id]["replaced_by"][0]
            return id
        if id in self.name_to_id:
            return self.name_to_id[id]
        if id in self.insdc_to_id:
            return self.insdc_to_id[id]
        raise ValueError(f"Unknown node: {id}")

    def node(self, id: str):
        nid = self.as_node_id(id)
        return self.graph.nodes[nid]

    def name_mapping(self):
        mapping = {}
        for so_id, node in self.nodes():
            mapping[so_id] = so_id
            name = node.get("name", None)
            if name:
                mapping[name] = so_id
            for insdc_name in insdc_synonyms(node):
                if insdc_name not in mapping:
                    mapping[insdc_name] = so_id
        return mapping

    def rna_type_tree(self, child, parents):
        if child in ALTERNATES:
            return ALTERNATES[child]

        print(child)
        for parent in parents:
            print(parent)
            if child == parent:
                break
            paths = nx.all_simple_paths(
                self.graph, source=child, target=parent)
            paths = list(paths)
            print(paths)
            if not paths:
                continue

            if len(paths) > 1:
                LOGGER.warn("Too many paths currently in %s", paths)

            tree = []
            for node_id in paths[0]:
                if node_id in SKIPPED_TERMS:
                    continue
                node = self.graph.nodes[node_id]
                tree.insert(0, (node_id, node["name"]))
            return tree

        LOGGER.error(
            "Assumes all SO terms are one of %s, %s is not", ", ".join(
                parents), child
        )
        return [(child, self.graph.nodes[child]["name"])]


@lru_cache()
def load_ontology(filename) -> SoOntology:
    return SoOntology.from_file(filename)


@lru_cache()
def rna_type_tree(ontology: SoOntology, child: str):
    nid = ontology.as_node_id(child)
    node = ontology.node(nid)
    if not node.get("so_term_tree", None):
        parents = [ontology.name_to_id[n] for n in BASE_SO_TERMS]
        node["so_term_tree"] = ontology.rna_type_tree(nid, parents)

    return node["so_term_tree"]


def insdc_synonyms(node):
    pattern = re.compile(r"INSDC_(?:feature|qualifier):(\w+)")
    for synonym in node.get("synonym", []):
        if "EXACT" not in synonym:
            continue
        match = re.search(pattern, synonym)
        if not match:
            continue
        yield match.group(1)


def name_index(ontology, filename) -> SqliteDict:
    mapping = SqliteDict(filename)
    mapping.update(ontology.name_mapping())
    mapping.commit()
    return mapping
