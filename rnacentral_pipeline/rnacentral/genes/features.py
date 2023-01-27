# -*- coding: utf-8 -*-

"""
Copyright [2009-${2022}] EMBL-European Bioinformatics Institute
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

import networkx as nx


def tree_distance(graph: nx.MultiDiGraph, node1: str, node2: str) -> int:
    r"""
    Calculate the distance through the graph from one node to another

        |            |
        |           / \
       / \         |   \
      /   \       / \   \
     A     B     C   D   E

    A-B distance should be 4 (2 up, 2 down)
    C-D distance should be 2 (1 up, 1 down)
    C-E distance should be 6 (3 up, 3 down)
    """

    path_1 = nx.descendants(graph, node1)
    path_2 = nx.descendants(graph, node2)
    uncommon_path = path_1.symmetric_difference(path_2) ## XOR operation

    return len(uncommon_path)


def exon_junction_jaccard_index(left: Region, right: Region) -> float:
    return 0.0
