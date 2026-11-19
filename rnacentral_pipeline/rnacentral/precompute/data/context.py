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

from pathlib import Path

import attr
import networkx as nx

from rnacentral_pipeline.databases.data import RnaType
from rnacentral_pipeline.databases.sequence_ontology import tree as so


@attr.s(frozen=True)
class Context:
    so_tree = attr.ib()

    @classmethod
    def from_directory(cls, path: Path) -> "Context":
        return cls(so_tree=so.load_ontology(so.REMOTE_ONTOLOGY))

    def so_term_for(self, name: str) -> RnaType:
        return RnaType.from_insdc_term(self.so_tree, name)

    def term_is_a(self, name: str, term: RnaType) -> bool:
        target = name
        if not name.startswith("SO:"):
            if name == "lnc_RNA":
                name = "lncRNA"
            target = self.so_tree.name_to_id[name]
        source = term.so_term.so_id
        if source == target:
            return True
        paths = nx.all_simple_paths(self.so_tree.graph, source=source, target=target)
        return bool(next(paths, None))
