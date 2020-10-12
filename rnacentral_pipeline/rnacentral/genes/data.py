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

import hashlib
import typing as ty

import attr
from intervaltree import Interval

from rnacentral_pipeline.databases.data.regions import Exon, SequenceRegion
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import (
    BedEntry, Strand)


@attr.s(auto_attribs=True, frozen=True)
class Extent:
    region_name: str
    assembly: str
    taxid: int
    chromosome: str
    strand: str
    start: int
    stop: int

    @classmethod
    def build(cls, raw):
        return cls(
            assembly=raw["assembly_id"],
            taxid=raw["taxid"],
            chromosome=raw["chromosome"],
            strand=raw["strand"],
            start=raw["start"],
            stop=raw["stop"],
        )

    def region_name(self, primary) -> str:
        return f'{primary}@{self.chromosome}/{self.start}-{self.stop}:{self.strand}'

    def merge(self, other: Extent):
        assert other.assembly == self.assembly
        assert other.taxid == self.taxid
        assert other.strand == self.strand
        assert other.chromosome = self.chromosome
        return attr.assoc(self, start=min(self.start, other.start),
                          stop=max(self.stop, other.stop))

    def as_region(self):
        return SequenceRegion(
            assembly_id=self.assembly,
            chromosome=self.chromosome,
            strand=Strand.build(self.strand),
            exons=Exon(start=self.start, stop=self.stop),
        )


@attr.s(auto_attribs=True, frozen=True)
class QaInfo:
    has_issue: bool

    @classmethod
    def build(cls, raw):
        return cls(has_issue=raw["has_issue"],)


@attr.s(auto_attribs=True, frozen=True)
class UnboundLocation:
    urs_taxid: str
    extent: Extent
    exons: ty.List[Exon]
    region_database_id: int
    insdc_rna_type: str
    so_rna_type: str
    qa: QaInfo
    providing_databases: ty.List[str]

    @classmethod
    def build(cls, raw):
        return cls(
            urs_taxid=raw["urs_taxid"],
            extent=Exten.build(raw),
            exons=[Exon.build(e) for e in raw["exons"]],
            region_database_id=raw["region_id"],
            insdc_rna_type=raw["insdc_rna_type"],
            so_rna_type=raw["so_rna_type"],
            qa=QaInfo.build(raw["qa"]),
            providing_databases=raw['providing_databases'],
        )

    @property
    def start(self):
        return self.extent.start

    @property
    def stop(self):
        return self.extent.stop

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self)


@attr.s(auto_attribs=True, frozen=True)
class LocusMember:
    info: UnboundLocation
    extent: Extent
    is_representative: bool

    @classmethod
    def from_location(cls, location):
        return cls(
            info=location,
            extent=location.extent,
            is_representative=False,
        )

    def as_bed(self):
        region = self.extent.as_region()
        region = attr.assoc(region, exons=self.exons)
        return BedEntry(
            rna_id=self.info.urs_taxid,
            rna_type=self.info.so_rna_type,
            databases=self.info.providing_databases,
            region=region,
        )


@attr.s(auto_attribs=True, frozen=True)
class Locus:
    extent: Extent
    members: ty.List[LocusMember]

    @classmethod
    def singleton(cls, location):
        return cls(
            extent=location.extent,
            members=[LocusMember.from_location(location)],
        )

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self)

    def merge(self, other: ty.List[Locus]) -> Locus:
        extent = self.extent
        mreged_members = set(self.members)
        for locus in other:
            extent.merge(locus.extent)
            members.update(locus.members)

        members = []
        merged_members = sorted(merged_members, key=op.attrgetter('extent'))
        for member in merged_members:
            members.append(attr.assoc(member, is_representative=False))

        return Locus(
            extent=extent,
            members=members
        )

    def id_hash(self):
        ids = {m.info.urs_taxid for m in self.members}
        ids = sorted(ids)
        return hashlib.sha256(''.join(ids))

    def locus_name(self):
        return self.extent.region_name(self.id_hash())

    def writeable(self):
        for member in self.members:
            yield [
                self.extent.taxid
                self.extent.assembly_id,
                self.locus_name,
                self.extent.chromosome,
                self.extent.strand,
                self.extent.start,
                self.extent.stop,
                member.urs_taxid,
                member.region_database_id,
                member.is_representative,
            ]

    def as_bed(self, include_gene=True, include_representaive=True,
               include_members=False) -> ty.List[BedEntry]:
        if include_gene:
            yield BedEntry(
                rna_id=self.id_hash(),
                rna_type='gene',
                databases='RNAcentral',
                region=self.extent.as_region(),
            )
        for member in self.members:
            if include_representative:
                if member.is_representative:
                    yield member.as_bed()
            elif include_members:
                yield member.as_bed()
