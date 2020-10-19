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

from . import data

REP_DBS = {
    "pdbe",
    "refseq",
    "ensembl",
    "hgnc",
}


def should_reject(location: data.UnboundLocation) -> bool:
    if any(d.lower() in REP_DBS for d in location.databases):
        return False
    return location.qa.has_issue


def intervals(locations: ty.Iterable[data.UnboundLocation]) -> IntervalTree:
    state = data.State()
    for location in locations:
        if should_reject(location):
            state.reject(location)
            continue
        locus = data.Locus.singleton(location)
        current = state.overlaps(location)
        if not current:
            state.add(locus)
        elif len(current) > 1:
            other = [i.data for i in current]
            current = locus.merge(other)
            state.add(current)
    return state


def mark_representative(members) -> ty.List[data.LocusMember]:
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


def build(locations) -> data.Finalized:
    state = intervals(locations)
    locuses = []
    for interval in state.tree:
        locus = interval.data
        updated = mark_representative(locus.members)
        locus = attr.assoc(locus, members=updated)
        locuses.append(locus)
    return data.Finalized(locuses, state.rejected)
