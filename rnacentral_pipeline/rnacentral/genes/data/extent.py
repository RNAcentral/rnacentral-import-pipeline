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

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data.regions import (CoordinateSystem, Exon,
                                                        SequenceRegion, Strand)


@attr.s(frozen=True, slots=True)
class Extent:
    region_name = attr.ib(validator=is_a(str))
    assembly = attr.ib(validator=is_a(str))
    taxid = attr.ib(validator=is_a(int))
    chromosome = attr.ib(validator=is_a(str))
    strand = attr.ib(validator=is_a(int))
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, raw):
        assert raw['region_start'] < raw['region_stop']
        return cls(
            assembly=raw["assembly_id"],
            taxid=raw["taxid"],
            chromosome=raw["chromosome"],
            strand=raw["strand"],
            start=raw["region_start"],
            stop=raw["region_stop"],
        )

    def string_strand(self):
        return Strand.build(self.strand).display_string()

    def region_name(self, primary) -> str:
        return f"{primary}@{self.chromosome}/{self.start}-{self.stop}:{self.strand}"

    def merge(self, other: "Extent"):
        assert other.assembly == self.assembly
        assert other.taxid == self.taxid
        assert other.strand == self.strand
        assert other.chromosome == self.chromosome
        updated = attr.assoc(
            self, start=min(self.start, other.start), stop=max(self.stop, other.stop)
        )
        assert updated.start < updated.stop
        return updated

    def as_region(self):
        return SequenceRegion(
            assembly_id=self.assembly,
            chromosome=self.chromosome,
            strand=Strand.build(self.strand),
            exons=[Exon(start=self.start, stop=self.stop)],
            coordinate_system=CoordinateSystem.zero_based(),
        )
