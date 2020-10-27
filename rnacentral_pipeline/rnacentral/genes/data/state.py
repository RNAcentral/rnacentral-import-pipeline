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
import logging
import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional
from intervaltree import IntervalTree

from rnacentral_pipeline.rnacentral.genes.data import (Cluster, ClusteringKey,
                                                       ClusterMember,
                                                       LocationInfo)

LOGGER = logging.getLogger(__name__)


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
    key = attr.ib(validator=is_a(ClusteringKey))
    clusters: ty.List[Cluster] = attr.ib(validator=is_a(list))
    rejected: ty.List[LocationInfo] = attr.ib(validator=is_a(list), factory=list)
    ignored: ty.List[LocationInfo] = attr.ib(validator=is_a(list), factory=list)

    def data_types(self):
        return [
            (DataType.clustered, self.clusters),
            (DataType.rejected, self.rejected),
            (DataType.ignored, self.ignored),
        ]


@attr.s(frozen=True, slots=True)
class LocationStatus:
    location: LocationInfo = attr.ib(validator=is_a(LocationInfo))
    status: ty.Optional[DataType] = attr.ib(validator=optional(is_a(DataType)))
    data: ty.Any = attr.ib()


@attr.s()
class State:
    key = attr.ib(validator=is_a(ClusteringKey))
    _tree = attr.ib(validator=is_a(IntervalTree), factory=IntervalTree)
    _locations: ty.Dict[int, LocationStatus] = attr.ib(
        validator=is_a(dict), factory=dict
    )
    _clusters: ty.Dict[int, Cluster] = attr.ib(validator=is_a(dict), factory=dict)

    def add_location(self, location: LocationInfo):
        LOGGER.debug("Adding location: %s", location.id)
        if location.id in self._locations:
            raise ValueError(f"Already seen location: {location}")
        self._locations[location.id] = LocationStatus(location, None, None)

    def overlaps(self, location: LocationInfo):
        intervals = self._tree.overlap(location.as_interval())
        LOGGER.debug(
            "Getting all overlaps for %s, found %i", location.id, len(intervals)
        )
        return [self._clusters[interval.data] for interval in intervals]

    def reject_location(self, location: LocationInfo):
        LOGGER.debug("Rejecting location %s", location.id)
        self.__remove_location_from_clusters__(location, (DataType.rejected, None))

    def ignore_location(self, location: LocationInfo):
        LOGGER.debug("Ignoring location %s", location.id)
        self.__remove_location_from_clusters__(location, (DataType.ignored, None))

    def add_to_cluster(self, location: LocationInfo, cluster_id: int):
        LOGGER.debug("Adding location %i to cluster %i", location.id, cluster_id)
        if cluster_id not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster_id}")

        cluster = self._clusters[cluster_id]
        del self._clusters[cluster_id]
        if cluster.as_interval() not in self._tree:
            raise ValueError(f"Cluster {cluster} is not indexed")
        self._tree.remove(cluster.as_interval())
        new_cluster = cluster.add_location(location)
        for member in new_cluster.members:
            location = member.location
            assert location.id in self._locations, "Unknown location {location}"
            del self._locations[location.id]
            self._locations[location.id] = LocationStatus(
                location, DataType.clustered, new_cluster.id
            )
        self._clusters[new_cluster.id] = new_cluster
        self._tree.add(new_cluster.as_interval())

    def merge_clusters(self, clusters: ty.List[Cluster]) -> int:
        assert len(clusters) > 1
        LOGGER.debug("Merging clusters: %s", [c.id for c in clusters])
        new_cluster = clusters[0]
        for cluster in clusters[1:]:
            assert cluster.id != new_cluster.id
            if cluster.id not in self._clusters:
                raise ValueError(f"Unknown cluster {cluster}")
            if cluster.as_interval() not in self._tree:
                raise ValueError(f"Unindexed cluster {cluster}")
            new_cluster = new_cluster.merge(cluster)
            del self._clusters[cluster.id]
            self._tree.remove(cluster.as_interval())

        self._clusters[new_cluster.id] = new_cluster
        self._tree.add(new_cluster.as_interval())
        for member in new_cluster.members:
            location = member.location
            if location.id not in self._locations:
                raise ValueError(f"Unknown location {location}")
            del self._locations[location.id]
            self._locations[location.id] = LocationStatus(
                location, DataType.clustered, new_cluster.id
            )
        assert new_cluster.id in self._clusters
        assert new_cluster.as_interval() in self._tree
        return new_cluster.id

    def add_singleton_cluster(self, location: LocationInfo):
        LOGGER.debug("Adding singleton cluster of %s", location.id)
        if location.id not in self._locations:
            raise ValueError(f"Unknown location: {location}")

        info = self._locations[location.id]
        if info.status is None or info.status in {
            DataType.rejected,
            DataType.ignored,
        }:
            cluster = Cluster.from_locations([location])
            LOGGER.info("Building singleton cluster %i", cluster.id)
            self._clusters[cluster.id] = cluster
            del self._locations[location.id]
            self._locations[location.id] = LocationStatus(
                location, DataType.clustered, cluster.id
            )
            self._tree.add(cluster.as_interval())
            assert cluster.id in self._clusters
            assert cluster.as_interval() in self._tree
        elif info.status == DataType.clustered:
            raise ValueError(f"Location {location} has already been clustered")
        else:
            raise ValueError(f"Unknown info {info} for {location}")

    def ignore_cluster(self, cluster_id: int):
        LOGGER.debug("Ignoring all locations in %s", cluster_id)
        if cluster_id not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster_id}")

        cluster = self._clusters[cluster_id]
        del self._clusters[cluster_id]
        if cluster.as_interval() not in self._tree:
            raise ValueError(f"Cluster {cluster} not indexed")

        self._tree.remove(cluster.as_interval())
        for member in cluster.members:
            location = member.location
            LOGGER.debug("Ignoring location %s", location.id)
            if location.id not in self._locations:
                raise ValueError(f"Somehow location {location} is unknown")
            del self._locations[location.id]
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
            key=self.key,
            clusters=list(self._clusters.values()),
            rejected=rejected,
            ignored=ignored,
        )

    def validate(self):
        assert len(self._clusters) <= len(self._locations)
        assert len(self._tree) == len(self._clusters)
        for cluster_id, cluster in self._clusters.items():
            assert cluster.as_interval() in self._tree

        for interval in self._tree:
            assert interval.data in self._clusters

        for location_id, info in self._locations.items():
            if info.data is None:
                LOGGER.debug("Skipping unprocessed location %s", info)
                continue
            if info.data not in self._clusters:
                raise ValueError("Location %s is in missing cluster", info)

    def __remove_location_from_clusters__(self, location: LocationInfo, new_status):
        if location.id not in self._locations:
            raise ValueError(f"Unknown location: {location}")

        info = self._locations[location.id]
        status = info.status
        if status is None or status in {DataType.ignored, DataType.rejected}:
            del self._locations[location.id]
            self._locations[location.id] = LocationStatus(location, *new_status)
        elif status == DataType.clustered:
            cluster = self._clusters[info.data]
            self._tree.remove(cluster.as_interval())
            new_cluster = cluster.remove_location(location)
            if new_cluster is None:
                del self._clusters[info.data]
                del self._locations[location.id]
                self._locations[location.id] = LocationStatus(location, *new_status)
            else:
                assert new_cluster.id == cluster.id
                del self._clusters[cluster.id]
                del self._locations[location.id]
                self._clusters[cluster.id] = cluster
                self._locations[location.id] = LocationStatus(location, *new_status)
                self._tree.add(new_cluster.as_interval())
        else:
            raise ValueError(f"Unknown status {status} for {location}")
