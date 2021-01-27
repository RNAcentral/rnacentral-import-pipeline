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
import tempfile
import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from intervaltree import Interval

from rnacentral_pipeline.rnacentral.repeats import ranges


@attr.s()
class RepeatTree:
    """
    This represents repeats across several assemblies.
    """

    store: ty.Dict[str, ranges.Info] = attr.ib(validator=is_a(dict), factory=dict)

    @classmethod
    def from_file(cls, path: Path) -> "RepeatTree":
        """
        Load the tree from the given file. This assumes that the format is the
        one used by `dump`.
        """
        with path.open("rb") as raw:
            tree = pickle.load(raw)
            tree.validate()
            return tree

    @classmethod
    def from_directory(cls, path: Path) -> "RepeatTree":
        return cls.from_file(path / "info.pickle")

    def add_info(self, info: ranges.Info):
        self.store[info.assembly_id] = info

    def has_assembly(self, assembly_id) -> bool:
        return assembly_id in self.store

    def assembly(self, assembly_id: str) -> ranges.RepeatRanges:
        return ranges.RepeatRanges.from_info(self.store[assembly_id])

    def validate(self):
        for info in self.store.values():
            info.validate()

    def dump(self, output: Path):
        """
        Write the tree to a file in a format suitable for RepeatTree.load.
        """
        assert output.is_dir(), f"{output} must be a directory"
        path = output / "info.pickle"
        with path.open("wb") as out:
            pickle.dump(attr.asdict(self), out)


def from_directories(paths: ty.List[Path], output: Path):
    tree = RepeatTree()
    for path in paths:
        assert path.is_dir()
        info = ranges.Info.from_directory(path)
        info.validate()
        tree.add_info(info)
    tree.validate()
    tree.dump(output)
