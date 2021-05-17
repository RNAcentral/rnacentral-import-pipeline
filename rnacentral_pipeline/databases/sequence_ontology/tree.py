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

import re
import logging
from functools import lru_cache

from sqlitedict import SqliteDict
import networkx as nx
import obonet


REMOTE_ONTOLOGY = 'https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so-simple.obo'

LOGGER = logging.getLogger(__name__)


@lru_cache()
def load_ontology(filename):
    ont = obonet.read_obo(filename)
    ont.id_to_name = {id_: data.get('name') for id_, data in ont.nodes(data=True)}
    ont.name_to_id = {data['name']: id_ for id_, data in ont.nodes(data=True) if 'name' in data}
    ont.insdc_to_id = {}
    for id_, node in ont.nodes(data=True):
        for insdc in insdc_synonyms(node):
            ont.insdc_to_id[insdc] = id_
    return ont


def compute_rna_type_tree(ontology, child, parent):
    paths = nx.all_simple_paths(ontology, source=child, target=parent)
    paths = list(paths)
    if not paths:
        LOGGER.error("Assumes all SO terms are ncRNA's currently: %s -> %s",
                     child, parent)
        return [(child, ontology.nodes[child]['name'])]

    if len(paths) > 1:
        LOGGER.warn("Too many paths currently in %s", paths)

    tree = []
    for node_id in paths[0]:
        node = ontology.nodes[node_id]
        tree.insert(0, (node_id, node['name']))
    return tree


@lru_cache()
def rna_type_tree(ontology, child):
    nid = None
    node = None
    if ontology.has_node(child):
        nid = child
        node = ontology.nodes[child]
    elif child in ontology.name_to_id:
        nid = ontology.name_to_id[child]
        node = ontology.nodes[nid]
    else:
        raise ValueError("Unknown node: " + child)

    base_node = ontology.name_to_id['ncRNA']
    if not node.get('so_term_tree', None):
        node['so_term_tree'] = compute_rna_type_tree(ontology, nid, base_node)
    return node['so_term_tree']


def insdc_synonyms(node):
    pattern = re.compile(r'INSDC_(?:feature|qualifier):(\w+)')
    for synonym in node.get('synonym', []):
        if 'EXACT' not in synonym:
            continue
        match = re.search(pattern, synonym)
        if not match:
            continue
        yield match.group(1)


def name_index(ontology, filename) -> SqliteDict:
    mapping = SqliteDict(filename)
    for so_id, node in ontology.nodes(data=True):
        mapping[so_id] = so_id
        name = node.get('name', None)
        if name:
            mapping[name] = so_id
        for insdc_name in insdc_synonyms(node):
            if insdc_name not in mapping:
                mapping[insdc_name] = so_id

    mapping.commit()
    return mapping
