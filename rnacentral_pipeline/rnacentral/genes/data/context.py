from __future__ import annotations

# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import collections as coll
import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from intervaltree import IntervalTree

from rnacentral_pipeline import psql
from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree

from .location import Count, LocationInfo


def load_counts(handle: ty.IO) -> ty.Dict[str, Count]:
    counts = {}
    for entry in psql.json_handler(handle):
        if entry["urs_taxid"] in counts:
            raise ValueError(f"Counts data contains duplicates: {entry['urs_taxid']}")
        counts[entry["urs_taxid"]] = Count.from_json(entry)
    return counts


def load_pseudogenes(path: Path) -> ty.DefaultDict[str, IntervalTree]:
    trees = coll.defaultdict(IntervalTree)
    with open(path, "r") as handle:
        for entry in psql.json_handler(handle):
            tree = trees[entry["chromosome"]]
            for exon in entry["exons"]:
                tree.addi(exon["exon_start"], exon["exon_stop"] + 1, entry["rna_id"])
            trees[entry["chromosome"]] = tree

    return trees


def load_repetitive(handle: ty.IO) -> ty.DefaultDict[str, IntervalTree]:
    trees = coll.defaultdict(IntervalTree)
    return trees


@attr.s()
class Context:
    ontology = attr.ib(validator=is_a(so_tree.SoOntology))
    pseudogenes: ty.DefaultDict[str, IntervalTree] = attr.ib(
        validator=is_a(ty.DefaultDict)
    )
    repetitive: ty.DefaultDict[str, IntervalTree] = attr.ib(
        validator=is_a(ty.DefaultDict)
    )
    counts: ty.Dict[str, Count] = attr.ib(validator=is_a(dict))
    max_rfam_shift = attr.ib(validator=is_a(int), default=10)

    @classmethod
    def from_files(cls, genes: Path, repetitive: ty.IO, counts: ty.IO) -> Context:
        ontology = so_tree.load_ontology(so_tree.REMOTE_ONTOLOGY)
        return cls(
            ontology=ontology,
            pseudogenes=load_pseudogenes(genes),
            repetitive=load_repetitive(repetitive),
            counts=load_counts(counts),
        )

    def count_for(self, urs_taxid: str) -> Count:
        if urs_taxid not in self.counts:
            raise ValueError(f"Missing counts for {urs_taxid}")
        return self.counts[urs_taxid]

    def overlaps_pseudogene(self, location: LocationInfo) -> bool:
        tree = self.pseudogenes[location.extent.chromosome]
        return bool(tree.overlaps(location.as_interval()))

    def overlaps_repetitive(self, location: LocationInfo) -> bool:
        tree = self.repetitive[location.extent.chromosome]
        return bool(tree.overlaps(location.as_interval()))
