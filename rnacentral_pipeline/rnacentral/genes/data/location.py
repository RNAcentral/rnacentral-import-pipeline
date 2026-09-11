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

import re
import typing as ty
from collections import OrderedDict

import attr
from attr.validators import instance_of as is_a
from gffutils import Feature
from intervaltree import Interval

from rnacentral_pipeline.databases.data.regions import Exon
from rnacentral_pipeline.databases.data.databases import Database
from rnacentral_pipeline.rnacentral.ftp_export.coordinates.bed import BedEntry

from .extent import Extent
from .rna_type import RnaType


@attr.s(auto_attribs=True, frozen=True, slots=True)
class QaInfo:
    has_issue: bool
    incomplete_sequence: bool
    possible_contamination: bool
    missing_rfam_match: bool

    @classmethod
    def build(cls, raw):
        return cls(
            has_issue=raw["has_issue"],
            incomplete_sequence=raw["incomplete_sequence"],
            possible_contamination=raw["possible_contamination"],
            missing_rfam_match=raw["missing_rfam_match"],
        )

    def as_tuple(self):
        return (
            self.has_issue,
            self.incomplete_sequence,
            self.possible_contamination,
            self.missing_rfam_match,
        )


@attr.s(frozen=True, slots=True, hash=True)
class Count:
    mapped_count = attr.ib(validator=is_a(int))
    given_count = attr.ib(validator=is_a(int))
    total_count = attr.ib(validator=is_a(int))

    @classmethod
    def from_json(cls, raw: ty.Dict[str, int]):
        return cls(
            mapped_count=raw["mapped_count"],
            given_count=raw["given_count"],
            total_count=raw["total_count"],
        )

    def __attrs_post_init__(self):
        if self.total_count <= 0:
            raise ValueError("Must have >0 total counts")
        if self.mapped_count + self.given_count != self.total_count:
            raise ValueError("Counts do not sum")


@attr.s(frozen=True, slots=True, hash=True)
class LocationInfo:
    id: int = attr.ib(validator=is_a(int))
    urs_taxid: str = attr.ib(validator=is_a(str))
    extent: Extent = attr.ib(validator=is_a(Extent))
    exons: ty.Tuple[Exon] = attr.ib(validator=is_a(tuple))
    region_name: str = attr.ib(validator=is_a(str))
    rna_type: RnaType = attr.ib(validator=is_a(RnaType))
    qa: QaInfo = attr.ib(validator=is_a(QaInfo))
    providing_databases: ty.Tuple[Database] = attr.ib(validator=is_a(tuple))
    databases: ty.Tuple[Database] = attr.ib(validator=is_a(tuple))
    counts: Count = attr.ib(validator=is_a(Count))

    @classmethod
    def build(cls, context, raw: ty.Dict[str, ty.Any]):
        rna_type = RnaType.build(
            raw["insdc_rna_type"], raw["so_rna_type"], context.ontology
        )
        return cls(
            id=raw["region_id"],
            urs_taxid=raw["urs_taxid"],
            extent=Extent.build(raw),
            exons=tuple(sorted(Exon.from_dict(e) for e in raw["exons"])),
            region_name=raw["region_name"],
            rna_type=rna_type,
            qa=QaInfo.build(raw["qa"][0]),
            providing_databases=tuple(
                Database.build(d) for d in raw["providing_databases"]
            ),
            databases=tuple(Database.build(d) for d in raw["databases"]),
            counts=context.count_for(raw["urs_taxid"]),
        )

    @property
    def start(self) -> int:
        return self.extent.start

    @property
    def stop(self) -> int:
        return self.extent.stop

    @property
    def name(self) -> str:
        return self.region_name

    @property
    def short_urs_suffix(self) -> str:
        urs, _taxid = self.urs_taxid.split("_")
        return re.sub("URS0*", "", urs)

    def has_introns(self):
        return len(self.exons) > 1

    def is_rfam_only(self):
        inferred = {Database.genecards, Database.malacards, Database.pirbase}
        dbs = set(self.databases) - inferred
        return dbs == set([Database.rfam])

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self.id)

    def as_bed(self) -> ty.Iterable[BedEntry]:
        region = self.extent.as_region()
        region = attr.assoc(region, exons=self.exons)
        yield BedEntry(
            rna_id=self.urs_taxid,
            rna_type=self.rna_type.so_term,
            databases=",".join(sorted(d.pretty() for d in self.providing_databases)),
            region=region,
        )

    def as_features(self) -> ty.Iterable[Feature]:
        transcript_name = self.urs_taxid
        yield Feature(
            seqid=self.extent.chromosome,
            source="RNAcentral",
            featuretype="transcript",
            start=self.extent.start,
            end=self.extent.stop,
            strand=self.extent.string_strand(),
            frame=".",
            attributes=OrderedDict(
                [
                    ("Name", [transcript_name]),
                    ("type", [self.rna_type.so_term]),
                ]
            ),
        )

        for index, exon in enumerate(self.exons):
            exon_name = transcript_name + f":ncRNA_exon{index + 1}"
            yield Feature(
                seqid=self.extent.chromosome,
                source="RNAcentral",
                featuretype="noncoding_exon",
                start=exon.start,
                end=exon.stop,
                strand=self.extent.string_strand(),
                frame=".",
                attributes=OrderedDict(
                    [
                        ("Name", [exon_name]),
                        ("Parent", [transcript_name]),
                        ("type", [self.rna_type.so_term]),
                    ]
                ),
            )

    def as_writeable(self, status=None) -> ty.Iterable[ty.List[str]]:
        yield [
            self.extent.assembly,
            self.id,
            self.urs_taxid,
            status,
        ]
