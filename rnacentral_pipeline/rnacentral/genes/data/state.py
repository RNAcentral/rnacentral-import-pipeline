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

import enum
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional
from intervaltree import IntervalTree

from . import Cluster, ClusterMember, LocationInfo


@enum.unique
class DataType(enum.Enum):
    ignored = enum.auto()
    rejected = enum.auto()
    clustered = enum.auto()

    @classmethod
    def all(cls):
        return [
            cls.clustered,
            cls.ignored,
            cls.rejected,
        ]


@attr.s(frozen=True)
class FinalizedState:
    chromosome = attr.ib(validator=is_a(str))
    clusters: ty.List[Cluster] = attr.ib(validator=is_a(list))
    rejected: ty.List[LocationInfo] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[LocationInfo] = attr.ib(validator=is_a(list), factory=list)

    def data_types(self):
        return [
            (DataType.clustered, self.clusters),
            (DataType.rejected, self.rejected),
            (DataType.ignored, self.ignored),
        ]


@attr.s(frozen=True)
class LocationStatus:
    location: LocationInfo = attr.ib(validator=is_a(LocationInfo))
    status: ty.Optional[DataType] = attr.ib(validator=optional(is_a(DataType)))
    data: ty.Any = attr.ib()


@attr.s()
class State:
    chromosome = attr.ib(validator=is_a(str))
    _tree = attr.ib(validator=is_a(IntervalTree), factory=IntervalTree)
    _locations: ty.Dict[int, LocationStatus] = attr.ib(
        validator=is_a(dict), factory=dict
    )
    _clusters: ty.Dict[int, Cluster] = attr.ib(validator=is_a(dict), factory=dict)

    def add_location(self, location: LocationInfo):
        if location.id in self._locations:
            raise ValueError(f"Already seen location: {location}")
        self._locations[location.id] = LocationStatus(location, None, None)

    def overlaps(self, location: LocationInfo):
        return self._tree.overlap(location.as_interval())

    def reject_location(self, location: LocationInfo):
        self.__remove_location_from_clusters__(location, (DataType.rejected, None))

    def ignore_location(self, location: LocationInfo):
        self.__remove_location_from_clusters__(location, (DataType.ignored, None))

    def merge_clusters(self, clusters: ty.List[int], additional=None):
        locations = additional or []
        tmp_cluster = Cluster.from_locations(locations)
        new_cluster = tmp_cluster.merge(clusters)
        self._clusters[new_cluster.id] = new_cluster
        self._tree.add(new_cluster.as_interval())

        for cluster in clusters:
            if cluster.id not in self._clusters:
                raise ValueError(f"Unknown cluster {cluster}")
            del self._clusters[cluster.id]
            self._tree.remove(cluster.as_interval())

        for member in new_cluster.members:
            location = member.location
            if location.id not in self._locations:
                raise ValueError(f"Unknown location {location}")
            self._locations[location.id] = LocationStatus(
                location, DataType.clustered, new_cluster.id
            )

    def add_singleton_cluster(self, location: LocationInfo):
        if location.id not in self._locations:
            raise ValueError(f"Unknown location: {location}")

        info = self._locations[location.id]
        if info.status is None or info.status in {
            DataType.rejected,
            DataType.ignored,
        }:
            cluster = Cluster.from_locations([location])
            self._clusters[cluster.id] = cluster
            self._locations[location.id] = LocationStatus(
                location, DataType.clustered, cluster.id
            )
            self._tree.add(cluster.as_interval())
        elif info.status == DataType.clustered:
            raise ValueError(f"Location {location} has already been clustered")
        else:
            raise ValueError(f"Unknown info {info} for {location}")

    def ignore_cluster(self, cluster_id: int):
        if cluster_id not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster_id}")

        cluster = self._clusters.pop(cluster_id)
        self._tree.remove(cluster.as_interval())
        for member in cluster.members:
            location = member.location
            if location.id not in self._locations:
                raise ValueError(f"Somehow location {location} is unknown")
            self._locations[location.id] = LocationStatus(
                location, DataType.ignored, None
            )

    def highlight_location(self, location_id: int):
        info = self._locations.get(location_id, None)
        if not info:
            raise ValueError(f"Unknown location id: {location_id}")

        if info.status != DataType.clustered:
            raise ValueError(f"Cannot highlight unclustered location: {info}")

        cluster = self._clusters.get(info.data, None)
        if not cluster:
            raise ValueError(f"Illegal State, missing cluster for {info}")

        new_cluster = cluster.highlight_location(info.location)
        assert new_cluster.id == cluster.id
        self._clusters[cluster.id] = new_cluster

    def members_of(self, cluster: int) -> ty.List[ClusterMember]:
        if cluster not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster}")
        return self._clusters[cluster].members

    def has_clusters(self) -> bool:
        return bool(self.clusters)

    def clusters(self) -> ty.List[int]:
        return list(self._clusters.keys())

    def finalize(self) -> FinalizedState:
        rejected: ty.List[LocationInfo] = []
        ignored: ty.List[LocationInfo] = []
        for (lid, data) in self._locations.items():
            if data.status is None:
                raise ValueError(f"Unclassified location {data.location}")
            elif data.status == DataType.clustered:
                continue
            elif data.status == DataType.rejected:
                rejected.append(data.location)
            elif data.status == DataType.ignored:
                ignored.append(data.location)
            else:
                raise ValueError(f"Unknown status state {data}")

        return FinalizedState(
            chromosome=self.chromosome,
            clusters=list(self._clusters.values()),
            rejected=rejected,
            ignored=ignored,
        )

    def __remove_location_from_clusters__(self, location, new_status):
        if location.id not in self._locations:
            raise ValueError(f"Unknown location: {location}")

        info = self._locations[location.id]
        status = info.status
        if status is None or status in {DataType.ignored, DataType.rejected}:
            self._locations[location.id] = LocationStatus(location, *new_status)
        elif status == DataType.clustered:
            cluster = self._clusters[info.data]
            self._tree.remove(cluster.as_interval())
            cluster.remove_location(location)
            if not cluster:
                del self._clusters[info.data]
            self._clusters[cluster.id] = cluster
            self._locations[location.id] = LocationStatus(location, *new_status)
            self._tree.add(cluster.as_interval())
        else:
            raise ValueError(f"Unknown status {status} for {location}")
