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

import pickle
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases.data import RnaType
from rnacentral_pipeline.rnacentral.repeats import tree
from rnacentral_pipeline.databases.sequence_ontology import tree as so


@attr.s(frozen=True)
class Context:
    so_tree = attr.ib(validator=is_a(SqliteDict))
    repeats = attr.ib(validator=is_a(tree.RepeatTree), factory=tree.RepeatTree)

    @classmethod
    def from_directory(cls, path: Path) -> "Context":
        repeats = tree.RepeatTree.from_directory(path / "repeat-tree")
        so_tree = so.load_ontology(so.REMOTE_ONTOLOGY)
        ctx = cls(repeats=repeats, so_tree=so_tree)
        ctx.validate()
        return ctx

    def validate(self):
        self.repeats.validate()

    def dump(self, path: Path):
        data = {'repeats': self.repeats}
        with path.open("wb") as out:
            pickle.dump(data, out)

    def so_term_for(self, name: str) -> RnaType:
        return RnaType.from_insdc_term(self.so_tree, name)
