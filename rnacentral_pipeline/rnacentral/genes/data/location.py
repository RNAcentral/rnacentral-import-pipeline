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
from collections import OrderedDict

import attr
from gffutils import Feature
from intervaltree import Interval

from rnacentral_pipeline.databases.data.regions import Exon
from rnacentral_pipeline.databases.sequence_ontology import tree as so_tree
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


@attr.s(auto_attribs=True, frozen=True, slots=True)
class UnboundLocation:
    urs_taxid: str
    extent: Extent
    exons: ty.Tuple[Exon]
    region_database_id: int
    rna_type: RnaType
    qa: QaInfo
    providing_databases: ty.Tuple[str]
    databases: ty.Tuple[str]

    @classmethod
    def build(cls, raw, ontology):
        rna_type = RnaType(
            insdc=raw["insdc_rna_type"],
            so_term=raw["so_rna_type"],
            ontology_terms=so_tree.rna_type_tree(ontology, raw["so_rna_type"]),
        )
        return cls(
            urs_taxid=raw["urs_taxid"],
            extent=Extent.build(raw),
            exons=tuple(sorted(Exon.from_dict(e) for e in raw["exons"])),
            region_database_id=raw["region_id"],
            rna_type=rna_type,
            qa=QaInfo.build(raw["qa"][0]),
            providing_databases=tuple(raw["providing_databases"]),
            databases=tuple(raw["databases"]),
        )

    @property
    def start(self):
        return self.extent.start

    @property
    def stop(self):
        return self.extent.stop

    def has_introns(self):
        return len(self.exons) > 1

    def is_rfam_only(self):
        return self.databases == ("Rfam",)

    def as_interval(self) -> Interval:
        return Interval(self.extent.start, self.extent.stop, self)

    def as_bed(self):
        region = self.extent.as_region()
        region = attr.assoc(region, exons=self.exons)
        return BedEntry(
            rna_id=self.urs_taxid,
            rna_type=self.so_rna_type,
            databases=",".join(self.providing_databases),
            region=region,
        )

    def as_features(self):
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
                [("Name", [transcript_name]), ("type", [self.so_rna_type]),]
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
                        ("type", [self.so_rna_type]),
                    ]
                ),
            )

        def as_writeable(self):
            return [
                self.extent.taxid,
                self.extent.assembly_id,
                self.locus_name,
                self.extent.chromosome,
                self.extent.strand,
                self.extent.start,
                self.extent.stop,
                self.urs_taxid,
                self.region_database_id,
                self.is_representative,
            ]
