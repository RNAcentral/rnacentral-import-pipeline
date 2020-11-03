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
import typing as ty
from pathlib import Path

import attr
import more_itertools as more
from attr.validators import instance_of as is_a
from Bio import SeqIO
from intervaltree import Interval, IntervalTree


@attr.s(frozen=True, slots=True)
class RepeatRanges:
    """
    This represents all repeat ranges in an assembly
    """

    assembly = attr.ib(validator=is_a(str))
    _chromosomes: ty.Dict[str, IntervalTree] = attr.ib(
        validator=is_a(dict), factory=dict
    )

    @classmethod
    def load(cls, path: Path) -> "RepeatRanges":
        """
        Load the RepeatRanges from the given path.
        """
        with path.open("rb") as raw:
            return pickle.load(raw)

    def insert(self, chromosome: str, start: int, stop: int):
        """
        Mark a given start/stop interval in a chromosome as a repetitive
        region.
        """

        if chromosome not in self._chromosomes:
            self._chromosomes[chromosome] = IntervalTree()
        self._chromosomes[chromosome].add(Interval(start, stop))

    def overlaps(self, chromosome: str, start: int, stop: int) -> ty.Set[Interval]:
        """
        Fetch all intervals that overlap the given start/stop in the given
        chromosome.
        """

        if chromosome not in self._chromosomes:
            raise ValueError(f"Unknown chromosome {chromosome}")
        tree = self._chromosomes[chromosome]
        return tree.overlaps(Interval(start, stop))

    def dump(self, output: Path):
        """
        Dump this object to the given path as a pickle file.
        """

        with output.open("wb") as out:
            pickle.dump(self, out)


def from_ensembl_fasta(assembly: str, path: Path) -> RepeatRanges:
    """
    Parse a fasta file from Ensembl that contains soft masked repeat information
    and produce a RepeatRanges representing that file. This assumes that the
    id for each entry in the file is the chromosome name.
    """

    ranges = RepeatRanges(assembly=assembly)
    for record in SeqIO.parse(path, "fasta"):
        chromosome = record.id
        nts = enumerate(record.seq, start=1)
        masked = more.split_at(nts, lambda nt: nt[1].upper() == nt[1])
        compressed = ((m[0][0], m[-1][0]) for m in masked)
        for (start, stop) in compressed:
            ranges.insert(chromosome, start, stop)
    return ranges
