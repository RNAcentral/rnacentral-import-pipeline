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
import operator as op
import typing as ty
import uuid
from collections import OrderedDict

import attr
from attr.validators import instance_of as is_a
from gffutils import Feature
from intervaltree import Interval

from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import BedEntry

from .extent import Extent
from .location import LocationInfo


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


@attr.s(auto_attribs=True, frozen=True, slots=True, hash=True)
class ClusterMember:
    location: LocationInfo
    member_type: MemberType = MemberType.member

    @classmethod
    def from_location(cls, location: LocationInfo):
        return cls(location=location)

    @property
    def extent(self):
        return self.location.extent

    def as_interval(self) -> Interval:
        return Interval(self.location.start, self.location.stop, self)

    def overlaps(self, other: "ClusterMember") -> bool:
        return self.as_interval().overlaps(other.as_interval())

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


@attr.s(slots=True)
class Cluster:
    extent: Extent = attr.ib(validator=is_a(Extent))
    members: ty.List[ClusterMember] = attr.ib(validator=is_a(list), factory=list)
    id: int = attr.ib(validator=is_a(int), factory=next_id)

    @classmethod
    def from_locations(cls, locations: ty.List[LocationInfo]) -> "Cluster":
        if not locations:
            raise ValueError("Cannot build with empty members")

        members = [ClusterMember.from_location(l) for l in locations]
        if len(locations) == 1:
            return cls(extent=locations[0].extent, members=members)
        extent = locations[0].extent
        for location in locations[1:]:
            extent = extent.merge(location.extent)
            members.append(ClusterMember.from_location(location))
        return cls(extent=extent, members=members,)

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self)

    def highlight_location(self, location: LocationInfo) -> "Cluster":
        updated = []
        for member in self.members:
            if member.location == location:
                updated.append(attr.evolve(member, member_type=MemberType.highlighted))
            else:
                updated.append(member)

        return attr.evolve(self, members=updated)

    def merge(self, other: ty.List["Cluster"]) -> "Cluster":
        extent = self.extent
        merged_members = set(self.members)
        for cluster in other:
            extent.merge(cluster.extent)
            merged_members.update(cluster.members)

        members = []
        sorted_members = sorted(merged_members, key=op.attrgetter("extent"))
        for member in sorted_members:
            members.append(attr.assoc(member, member_type=MemberType.member))

        return Cluster(extent=extent, members=members)

    def id_hash(self):
        ids = {m.location.urs_taxid for m in self.members}
        ids = sorted(ids)
        ids = "".join(ids)
        hash = hashlib.sha256(ids.encode("ascii", "ignore"))
        return hash.hexdigest()

    @property
    def name(self):
        return self.cluster_name()

    def rna_types(self):
        types = {m.location.rna_type for m in self.members}
        return types

    def cluster_name(self):
        return self.extent.region_name(self.id_hash())

    def remove_member(self, location):
        pass

    def validate_members(self):
        if not self.members:
            raise EmptyCluster(self)
        if len(self.members) == 1:
            return True
        return True

    def locations(self) -> ty.List[LocationInfo]:
        return [m.location for m in self.members]

    def __writeable_members__(self, allowed_members) -> ty.List[ClusterMember]:
        return [m for m in self.members if m.member_type in allowed_members]

    def as_writeable(self, allowed_members={MemberType.highlighted}):
        members = self.__writeable_members__(allowed_members)
        for member in members:
            yield [
                self.extent.taxid,
                self.extent.assembly,
                self.cluster_name(),
                self.extent.chromosome,
                self.extent.strand,
                self.extent.start,
                self.extent.stop,
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
