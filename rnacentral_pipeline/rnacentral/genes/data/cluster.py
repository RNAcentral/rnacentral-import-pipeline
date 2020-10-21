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
from collections import OrderedDict

import attr
from intervaltree import Interval, IntervalTree
from gffutils import Feature

from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import BedEntry

from .extent import Extent
from .location import UnboundLocation


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


@enum.unique
class MemberType(enum.Enum):
    highlighted = enum.auto()
    member = enum.auto()


@attr.s(auto_attribs=True, frozen=True, slots=True, hash=True)
class ClusterMember:
    info: UnboundLocation
    extent: Extent
    member_type: MemberType = MemberType.member

    @classmethod
    def from_location(cls, location):
        return cls(info=location, extent=location.extent)

    def as_interval(self) -> Interval:
        return Interval(self.info.start, self.info.stop, self)

    def overlaps(cls, other: 'ClusterMember') -> bool:
        return self.as_interval().overlaps(other.as_interval())

    def as_bed(self) -> ty.Iterable[BedEntry]:
        yield self.info.as_bed()

    def as_features(self) -> ty.Iterable[Feature]:
        features = self.info.as_features()
        transcript = next(features)
        transcript.attributes["member_type"] = [str(self.member_type.name.lower())]
        yield transcript
        for feature in features:
            yield feature

    def as_writeable(self):
        return self.info.writeable()


@attr.s(auto_attribs=True, frozen=True, slots=True)
class Cluster:
    extent: Extent
    members: ty.List[ClusterMember]

    @classmethod
    def singleton(cls, location):
        return cls(
            extent=location.extent, members=[ClusterMember.from_location(location)],
        )

    @classmethod
    def from_locations(cls, members: ty.List[UnboundLocation]) -> 'Cluster':
        if not members:
            raise ValueError("Cannot build with empty members")
        if len(members) == 1:
            return cls.singleton(members[0])
        cluster = cls.singleton(members[0])
        cluster = cluster.merge(members[1:])
        return cluster

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self)

    def merge(self, other: ty.List["Cluster"]) -> "Cluster":
        extent = self.extent
        merged_members = set(self.members)
        for cluster in other:
            extent.merge(cluster.extent)
            merged_members.update(cluster.members)

        members = []
        merged_members = sorted(merged_members, key=op.attrgetter("extent"))
        for member in merged_members:
            members.append(attr.assoc(member, is_representative=False))

        return Cluster(extent=extent, members=members)

    def id_hash(self):
        ids = {m.info.urs_taxid for m in self.members}
        ids = sorted(ids)
        ids = "".join(ids)
        hash = hashlib.sha256(ids.encode("ascii", "ignore"))
        return hash.hexdigest()

    def cluster_name(self):
        return self.extent.region_name(self.id_hash())

    def validate_members(self):
        if not self.members:
            raise EmptyCluster(self)
        if len(self.members) == 1:
            return True
        first = self.members[0]
        for member in self.members[1:]:
            if not member.overlaps(first):
                raise NonOverlappingMember(self)
        return True

    def __writeable_members__(self, allowed_members) -> ty.List[ClusterMember]:
        return [m for m in self.members if m.member_type in allowed_members]

    def as_writeable(self, allowed_members={MemberType.highlighted}):
        members = self.__writeable_members__(allowed_members)
        for member in members:
            yield [
                self.extent.taxid,
                self.extent.assembly_id,
                self.cluster_name,
                self.extent.chromosome,
                self.extent.strand,
                self.extent.start,
                self.extent.stop,
                member.urs_taxid,
                member.region_database_id,
                member.is_representative,
            ]

    def as_bed(
        self, include_gene=True, allowed_members={MemberType.highlighted}
    ) -> ty.List[BedEntry]:
        written_gene = not include_gene
        members = self.__writeable_members__(allowed_members)
        if members and include_gene:
            yield BedEntry(
                rna_id=self.id_hash(),
                rna_type="gene",
                databases="RNAcentral",
                region=self.extent.as_region(),
            )
        for member in members:
            yield member.as_bed()

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
                attributes=OrderedDict([("ID", [gene_id]),]),
            )
        for member in self.members:
            for feature in member.as_features():
                if feature.featuretype == "transcript":
                    feature.attributes["Parent"] = [gene_id]
                yield feature
