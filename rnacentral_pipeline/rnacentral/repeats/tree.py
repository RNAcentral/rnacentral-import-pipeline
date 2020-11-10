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

import shelve
import tempfile
import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from intervaltree import Interval

from rnacentral_pipeline.rnacentral.repeats import ranges


class RepeatTree:
    """
    This represents repeats across several assemblies.
    """

    def __init__(self, output: ty.Optional[Path]):
        if output:
            filename = str(output)
        else:
            filename = tempfile.mktemp()
        self.store = shelve.open(filename)

    @classmethod
    def load(cls, path: Path) -> "RepeatTree":
        """
        Load the tree from the given file. This assumes that the format is the
        one used by `dump`.
        """
        return cls(path)

    def has_assembly(self, assembly_id):
        return assembly_id in self.store

    def overlaps(
        self, assembly: str, chromosome: str, start: int, stop: int
    ) -> ty.Set[Interval]:
        """
        Find the overlaps for the given assembly/chromosome/start/stop.
        """

        if not self.has_assembly(assembly):
            raise ValueError(f"Unknown assembly {assembly}")

        return self.store[assembly].overlaps(chromosome, start, stop)

    def envelops(
        self, assembly: str, chromosome: str, start: int, stop: int
    ) -> ty.Set[Interval]:
        """
        Find all regions that enclose the given assembly/chromosome/start/stop.
        """

        if not self.has_assembly(assembly):
            raise ValueError(f"Unknown assembly {assembly}")

        return self.store[assembly].envelops(chromosome, start, stop)

    def add(self, value: ranges.RepeatRanges):
        """
        Add a given RepeatRange to this tree. This fails if the assembly is
        already stored in the tree.
        """

        if self.has_assembly(value.assembly):
            raise ValueError(f"Duplicate assmebly {value.assembly}")
        self.store[value.assembly] = value

    def dump(self):
        """
        Write the tree to a file in a format suitable for RepeatTree.load.
        """
        self.store.sync()
        self.close()


def from_ranges(output: Path, paths: ty.List[Path]) -> RepeatTree:
    """
    Build a repeat tree from a list of file handles that contain individual
    RepeatRanges.
    """

    tree = RepeatTree(output)
    for path in paths:
        loaded = ranges.RepeatRanges.load(path)
        tree.add(loaded)
        del loaded
    return tree
