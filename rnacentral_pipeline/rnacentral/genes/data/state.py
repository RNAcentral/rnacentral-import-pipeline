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
from intervaltree import Interval, IntervalTree

from . import Cluster, UnboundLocation


@attr.s(frozen=True, auto_attribs=True)
class StateUpdate:
    reject: ty.Set[UnboundLocation] = set()
    ignore: ty.Set[UnboundLocation] = set()
    new_clusters: ty.List[Cluster] = []
    old_clusters: ty.List[Cluster] = []

    @classmethod
    def no_op(cls):
        return cls(reject={}, ignore={}, new_clusters=[], old_clusters=[])

    def validate(self):
        for new_cluster in self.new_clusters:
            new_cluster.validate_members()


@attr.s(frozen=True)
class FinalizedState:
    chromosome = attr.ib(validator=is_a(str))
    clusters: ty.List[Cluster] = attr.ib(validator=is_a(list))
    rejected: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)

    def data_types(self):
        return [
            ("genes", self.clusters),
            ("rejected", self.rejected),
            ("ignored", self.ignored),
        ]


@attr.s()
class State:
    chromosome = attr.ib(validator=is_a(str))
    tree = attr.ib(validator=is_a(IntervalTree), factory=IntervalTree)
    rejected: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)

    def overlaps(self, location: UnboundLocation):
        return self.tree.overlap(location.start, location.stop)

    def reject(self, location: UnboundLocation):
        self.rejected.append(location)

    def ignore(self, location: UnboundLocation):
        self.ignored.append(location)

    def remove_interval(self, interval: Interval):
        self.tree.remove(interval)

    def update(self, update):
        self.rejected.extend(update.reject)
        self.ignored.extend(update.ignore)

        for old in update.old_clusters:
            tree.remove(old.as_interval())
        tree.update(u.as_interval for u in update.new_clusters)

    def has_clusters(self) -> bool:
        return bool(len(self.tree))

    def clusters(self) -> ty.Iterable[Cluster]:
        for interval in self.tree:
            yield interval.data

    def finalize(self) -> FinalizedState:
        return FinalizedState(
            chromosome=self.chromosome,
            clusters=list(self.tree),
            rejected=self.rejected,
            ignored=self.ignored,
        )
