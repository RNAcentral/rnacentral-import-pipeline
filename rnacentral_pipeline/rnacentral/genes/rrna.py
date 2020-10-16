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

import typing as ty

import attr
from attr.validators import instance_of as is_a
from intervaltree import IntervalTree

from .data import Locus, LocusMember, UnboundLocation

REP_DBS = {
    "pdbe",
    "refseq",
    "ensembl",
    "hgnc",
}


@attr.s()
class State:
    tree = attr.ib(validator=is_a(IntervalTree), factory=IntervalTree)
    rejected: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)

    def reject(self, location: UnboundLocation):
        self.rejected.append(location)

    def overlaps(self, location: UnboundLocation):
        return self.tree.overlap(location.start, location.stop)

    def add(self, locus: Locus):
        self.tree.add(locus.as_interval())


def should_reject(location: UnboundLocation) -> bool:
    if "PDBe" in location.databases:
        return False
    return location.qa.has_issue


def intervals(locations: ty.Iterable[UnboundLocation]) -> IntervalTree:
    state = State()
    for location in locations:
        if should_reject(location):
            state.reject(location)
            continue
        locus = Locus.singleton(location)
        current = state.overlaps(location)
        if not current:
            state.add(locus)
        elif len(current) > 1:
            other = [i.data for i in current]
            current = locus.merge(other)
            state.add(current)
    return state


def mark_representative(members) -> ty.List[LocusMember]:
    updates = []
    for member in members:
        qa = member.info.qa
        is_representative = False
        if any(d.lower() in REP_DBS for d in member.info.databases):
            is_representative = True
        elif not qa.has_issue and not member.info.has_introns():
            is_representative = True
        updates.append(attr.assoc(member, is_representative=is_representative))
    return updates


def build(locations):
    state = intervals(locations)
    for interval in state.tree:
        locus = interval.data
        updated = mark_representative(locus.members)
        locus = attr.assoc(locus, members=updated)
        yield locus
