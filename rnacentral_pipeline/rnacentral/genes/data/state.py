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


@attr.s(frozen=True)
class FinalizedState:
    chromosome = attr.ib(validator=is_a(str))
    clusters: ty.List[Cluster] = attr.ib(validator=is_a(list))
    rejected: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)


@attr.s()
class State:
    chromosome = attr.ib(validator=is_a(str))
    tree = attr.ib(validator=is_a(IntervalTree), factory=IntervalTree)
    rejected: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[UnboundLocation] = attr.ib(validator=is_a(list), factory=list)

    def reject(self, location: UnboundLocation):
        self.rejected.append(location)

    def reject_all(self, locations: UnboundLocation):
        self.rejected.extend(locations)

    def ignore(self, location: UnboundLocation):
        self.ignored.append(location)

    def ignore_all(self, locations: UnboundLocation):
        self.ignored.extend(locations)

    def overlaps(self, location: UnboundLocation):
        return self.tree.overlap(location.start, location.stop)

    def remove_interval(self, interval: Interval):
        self.tree.remove(interval)

    def add_clusters(self, clusters: ty.List[Cluster]):
        for cluster in clusters:
            self.tree.add(cluster.as_interval())

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
