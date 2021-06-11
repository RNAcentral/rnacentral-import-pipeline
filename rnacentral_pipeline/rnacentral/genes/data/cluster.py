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
import hashlib
import logging
import typing as ty
import uuid
from collections import OrderedDict

import attr
from attr.validators import instance_of as is_a
from gffutils import Feature
from intervaltree import Interval

from rnacentral_pipeline.databases.helpers.hashes import crc64
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import BedEntry

from .extent import Extent
from .location import LocationInfo
from .rna_type import RnaType

LOGGER = logging.getLogger(__name__)


class EmptyCluster(Exception):
    """
    Raised if there is an empty Cluster for some reason.
    """

    pass


class NonOverlappingMember(Exception):
    """
    Raised if a cluster contains members that do not overlap.
    """

    pass


def next_id():
    return uuid.uuid4().int


@enum.unique
class MemberType(enum.Enum):
    highlighted = enum.auto()
    member = enum.auto()


@attr.s()
class ClusteringKey:
    chromosome = attr.ib(validator=is_a(str))
    strand = attr.ib(validator=is_a(int))

    @classmethod
    def from_location(cls, location: LocationInfo):
        return cls(chromosome=location.extent.chromosome, strand=location.extent.strand)


@attr.s(slots=True)
class ClusterMember:
    location_id = attr.ib(validator=is_a(int))
    rna_type = attr.ib(validator=is_a(RnaType))
    extent = attr.ib(validator=is_a(Extent))
    member_type = attr.ib(validator=is_a(MemberType))

    @classmethod
    def from_location(
        cls, location: LocationInfo, member_type=MemberType.member
    ) -> "ClusterMember":
        return cls(
            location_id=location.id,
            rna_type=location.rna_type,
            extent=location.extent,
            member_type=member_type,
        )

    @classmethod
    def from_member(
        cls, member: "ClusterMember", member_type=MemberType.member
    ) -> "ClusterMember":
        return cls(
            location_id=member.location_id,
            rna_type=member.rna_type,
            extent=member.extent,
            member_type=member_type,
        )


@attr.s(slots=True)
class Cluster:
    extent: Extent = attr.ib(validator=is_a(Extent))
    _members: ty.Dict[int, ClusterMember] = attr.ib(validator=is_a(dict), factory=dict)
    id: int = attr.ib(validator=is_a(int), factory=next_id)

    @classmethod
    def from_locations(cls, locations: ty.List[LocationInfo]) -> "Cluster":
        if not locations:
            raise ValueError("Cannot build with empty members")

        members = {l.id: ClusterMember.from_location(l) for l in locations}
        if len(members) != len(locations):
            raise ValueError("Cannot build cluster with duplicate locations")

        if len(locations) == 1:
            return cls(extent=locations[0].extent, members=members)

        extent = locations[0].extent
        for location in locations[1:]:
            extent = extent.merge(location.extent)
        return cls(
            extent=extent,
            members=members,
        )

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self.id)

    def remove_location(self, location: LocationInfo):
        if location.id not in self._members:
            raise ValueError(f"Location {location} not part of cluster")

        if len(self._members) == 1:
            raise ValueError("Cannot create empty cluster")

        del self._members[location.id]
        first = next(iter(self._members.values()))
        self.extent = first.extent
        for member in self._members.values():
            self.extent = self.extent.merge(member.extent)

    def add_location(self, location: LocationInfo):
        if location.id in self._members:
            raise ValueError(f"Cannot add duplicate location {location}")

        self._members[location.id] = ClusterMember.from_location(location)
        self.extent = self.extent.merge(location.extent)

    def highlight_location(self, location: LocationInfo):
        if location.id not in self._members:
            raise ValueError(f"Unknown location {location}")

        self._members[location.id].member_type = MemberType.highlighted

    def merge(self, cluster: "Cluster"):
        for lid in self._members.keys():
            self._members[lid].member_type = MemberType.member

        for lid, member in cluster._members.items():
            if lid in self._members:
                raise ValueError(f"Illegal state, location {lid} in two clusters")
            self._members[lid] = ClusterMember.from_member(member)
            self.extent = self.extent.merge(member.extent)

    def rna_types(self) -> ty.Set[RnaType]:
        return {m.rna_type for m in self._members.values()}

    def location_ids(self) -> ty.List[int]:
        return list(self._members.keys())

    def __len__(self):
        return len(self._members)


@attr.s(frozen=True, slots=True)
class WriteableClusterMember:
    location = attr.ib(validator=is_a(LocationInfo))
    member_type = attr.ib(validator=is_a(MemberType))

    @classmethod
    def build(
        cls, member: ClusterMember, location: LocationInfo
    ) -> "WriteableClusterMember":
        return cls(location=location, member_type=member.member_type)

    def as_bed(self) -> ty.Iterable[BedEntry]:
        yield from self.location.as_bed()

    def as_features(self) -> ty.Iterable[Feature]:
        features = self.location.as_features()
        transcript = next(features)
        transcript.attributes["member_type"] = [str(self.member_type.name.lower())]
        yield transcript
        for feature in features:
            yield feature

    def as_writeable(self):
        return self.location.writeable()


@attr.s(slots=True, frozen=True)
class WriteableCluster:
    extent: Extent = attr.ib(validator=is_a(Extent))
    members: ty.Set[WriteableClusterMember] = attr.ib(validator=is_a(set), factory=set)

    @classmethod
    def build(
        cls, cluster: Cluster, mapping: ty.Dict[int, LocationInfo]
    ) -> "WriteableCluster":
        members = set()
        for lid, member in cluster._members.items():
            location = mapping[lid]
            members.add(WriteableClusterMember.build(member, location))
        return cls(extent=cluster.extent, members=members)

    def id_hash(self):
        ids = {m.location.urs_taxid for m in self.members}
        ids = sorted(ids)
        ids = "".join(ids)
        hash = hashlib.sha256(ids.encode("ascii", "ignore"))
        return hash.hexdigest()

    @property
    def name(self):
        return self.cluster_name()

    def cluster_name(self):
        return self.extent.region_name(self.id_hash())

    def cluster_hash(self):
        return crc64(self.cluster_name())

    def __writeable_members__(self, allowed_members) -> ty.List[WriteableClusterMember]:
        return [m for m in self.members if m.member_type in allowed_members]

    def as_writeable(self, allowed_members={MemberType.highlighted}):
        members = self.__writeable_members__(allowed_members)
        for member in members:
            yield [
                self.extent.assembly,
                self.cluster_name(),
                self.cluster_hash(),
                self.extent.chromosome,
                self.extent.strand,
                self.extent.start,
                self.extent.stop,
                len(self.members),
                member.location.urs_taxid,
                member.location.id,
                member.member_type.name,
            ]

    def as_bed(
        self, include_gene=True, allowed_members={MemberType.highlighted}
    ) -> ty.Iterable[BedEntry]:
        members = self.__writeable_members__(allowed_members)
        if members and include_gene:
            yield BedEntry(
                rna_id=self.id_hash(),
                rna_type="gene",
                databases="RNAcentral",
                region=self.extent.as_region(),
            )
        for member in members:
            yield from member.as_bed()

    def as_features(
        self, include_gene=True, allowed_members={MemberType.highlighted}
    ) -> ty.Iterable[Feature]:
        gene_id = self.id_hash()
        members = self.__writeable_members__(allowed_members)
        if members and include_gene:
            yield Feature(
                seqid=self.extent.chromosome,
                source="RNAcentral",
                featuretype="gene",
                start=self.extent.start,
                end=self.extent.stop,
                strand=self.extent.string_strand(),
                frame=".",
                attributes=OrderedDict([("ID", [gene_id])]),
            )
        for member in self.members:
            for feature in member.as_features():
                if feature.featuretype == "transcript":
                    feature.attributes["Parent"] = [gene_id]
                yield feature

    def __len__(self):
        return len(self.members)
