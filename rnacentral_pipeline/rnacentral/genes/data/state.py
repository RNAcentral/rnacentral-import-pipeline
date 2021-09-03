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

from rnacentral_pipeline.rnacentral.genes.data import (
    Cluster,
    ClusteringKey,
    LocationInfo,
    WriteableCluster,
)

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
    clusters: ty.List[WriteableCluster] = attr.ib(validator=is_a(list))
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
    method = attr.ib(validator=is_a(str))
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
        self.__remove_location_from_clusters__(location, DataType.rejected)
        assert len(self._clusters) == len(self._tree)

    def ignore_location(self, location: LocationInfo):
        LOGGER.debug("Ignoring location %s", location.id)
        self.__remove_location_from_clusters__(location, DataType.ignored)
        assert len(self._clusters) == len(self._tree)

    def add_to_cluster(self, location: LocationInfo, cluster_id: int):
        LOGGER.debug("Updating cluster %i with location %i", cluster_id, location.id)
        if cluster_id not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster_id}")

        if location.id not in self._locations:
            raise ValueError(f"Unknown location {location}")

        cluster = self._clusters[cluster_id]
        del self._clusters[cluster_id]
        interval = cluster.as_interval()
        if interval not in self._tree:
            raise ValueError(f"Cluster {cluster} is not indexed")
        self._tree.remove(interval)
        assert len(self._clusters) == len(self._tree)

        cluster.add_location(location)
        self._locations[location.id] = LocationStatus(
            location, DataType.clustered, cluster.id
        )
        self._clusters[cluster.id] = cluster
        self._tree.add(cluster.as_interval())
        if len(self._tree) != len(self._clusters):
            raise ValueError(
                f"Tree and cluster mismatch {len(self._tree)} vs {len(self._clusters)}"
            )

    def merge_clusters(self, clusters: ty.List[Cluster]) -> int:
        assert len(clusters) > 1
        new_cluster = clusters[0]
        if new_cluster.id not in self._clusters:
            raise ValueError(f"Unknown cluster {new_cluster}")
        if new_cluster.as_interval() not in self._tree:
            raise ValueError(f"Unindexed cluster {new_cluster}")
        del self._clusters[new_cluster.id]
        self._tree.remove(new_cluster.as_interval())

        for cluster in clusters[1:]:
            LOGGER.debug("Merging cluster {cluster.id} into {new_cluster.id}")
            if cluster.id not in self._clusters:
                raise ValueError(f"Unknown cluster {cluster}")
            if cluster.as_interval() not in self._tree:
                raise ValueError(f"Unindexed cluster {cluster}")
            new_cluster.merge(cluster)
            del self._clusters[cluster.id]
            self._tree.remove(cluster.as_interval())
            for lid in cluster._members.keys():
                if lid not in self._locations:
                    raise ValueError(f"Unknown location {lid} in cluster {cluster}")
                info = self._locations[lid]
                assert info.data == cluster.id
                del self._locations[lid]
                self._locations[info.location.id] = LocationStatus(
                    info.location, DataType.clustered, new_cluster.id
                )
        self._clusters[new_cluster.id] = new_cluster
        self._tree.add(new_cluster.as_interval())
        assert len(self._clusters) == len(self._tree)
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
        for lid, member in cluster._members.items():
            if lid not in self._locations:
                raise ValueError(f"Somehow location {lid} is unknown")
            info = self._locations[lid]
            del self._locations[lid]
            LOGGER.debug("Ignoring location %s", info.location.id)
            self._locations[lid] = LocationStatus(info.location, DataType.ignored, None)

    def highlight_location(self, location_id: int):
        info = self._locations.get(location_id, None)
        if not info:
            raise ValueError(f"Unknown location id: {location_id}")

        if info.status != DataType.clustered:
            raise ValueError(f"Cannot highlight unclustered location: {info}")

        cluster = self._clusters.get(info.data, None)
        if not cluster:
            raise ValueError(f"Illegal State, missing cluster for {info}")

        cluster.highlight_location(info.location)

    def members_of(self, cluster_id: int) -> ty.List[LocationInfo]:
        if cluster_id not in self._clusters:
            raise ValueError(f"Unknown cluster {cluster_id}")
        members = []
        cluster = self._clusters[cluster_id]
        for lid in cluster.location_ids():
            if lid not in self._locations:
                raise ValueError("Unknown location id {lid} in cluster {cluster}")
            info = self._locations[lid]
            if info.status != DataType.clustered or info.data != cluster_id:
                raise ValueError(f"Info {info} does not show clustered to {cluster_id}")
            members.append(info.location)
        return members

    def has_clusters(self) -> bool:
        return bool(self.clusters)

    def has_cluster(self, cluster_id: int) -> bool:
        return cluster_id in self._clusters

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

        clusters: ty.List[WriteableCluster] = []
        for cid, cluster in self._clusters.items():
            locations = {
                lid: self._locations[lid].location for lid in cluster.location_ids()
            }
            clusters.append(WriteableCluster.build(cluster, locations))

        return FinalizedState(
            key=self.key,
            clusters=clusters,
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

    def lengths(self):
        return (len(self._locations), len(self._tree), len(self._clusters))

    def __remove_location_from_clusters__(
        self, location: LocationInfo, new_status: DataType
    ):
        if location.id not in self._locations:
            raise ValueError(f"Unknown location: {location}")

        info = self._locations[location.id]
        del self._locations[location.id]

        status = info.status
        if status is None or status in {DataType.ignored, DataType.rejected}:
            self._locations[location.id] = LocationStatus(location, new_status, None)

        elif status == DataType.clustered:
            cluster = self._clusters[info.data]
            self._locations[location.id] = LocationStatus(location, new_status, None)
            if len(cluster) == 1:
                LOGGER.debug(
                    f"Cluster {info.data} has only one member, removing it completely"
                )
                del self._clusters[info.data]
                self._tree.remove(cluster.as_interval())
                assert info.data not in self._clusters
            else:
                assert location.id in cluster._members
                del self._clusters[info.data]
                self._tree.remove(cluster.as_interval())
                cluster.remove_location(location)
                self._tree.add(cluster.as_interval())
                self._clusters[cluster.id] = cluster
                assert location.id not in cluster._members
        else:
            raise ValueError(f"Unknown status {status} for {location}")
        assert len(self._tree) == len(self._clusters)
