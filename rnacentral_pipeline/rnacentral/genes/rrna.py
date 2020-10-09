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

import attr
from intervaltree import IntervalTree


def intervals(locations) -> IntervalTree:
    tree = IntervalTree()
    for location in locations:
        current = tree.overlap(location.start, location.stop)
        locus_interval = Locus.singleton(location)
        if current:
            other = [i.data for i in current]
            locus = locus_interval.merge(other)
        tree.add(locus.interval())
    return tree


def mark_representative(members):
    updates = []
    for member in members:
        qa = member.info.qa
        is_representative = not qa.has_issue
        updates.append(attr.assoc(member, is_representative=is_representative))
    return updates


def build(locations):
    tree = intervals(locations)
    for interval in tree:
        locus = interval.data
        updated = mark_representative(locus.members)
        locus = attr.assoc(locus, members=updated)
        yield locus
