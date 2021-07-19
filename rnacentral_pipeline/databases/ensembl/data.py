# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.data import SequenceRegion
from rnacentral_pipeline.databases.data.utils import as_so_term
from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.ensembl import helpers

LOGGER = logging.getLogger(__name__)


@enum.unique
class Division(enum.Enum):
    fungi = enum.auto()
    plants = enum.auto()
    protists = enum.auto()
    metazoa = enum.auto()
    vertebrates = enum.auto()

    @classmethod
    def from_name(cls, name):
        for value in cls:
            if value.name.lower() == name.lower():
                return value
        raise ValueError("Unknown Division %s" % name)

    @classmethod
    def names(cls):
        return [x.name for x in cls]

    @property
    def division_name(self):
        name = self.name[0].upper() + self.name[1:]
        return f"Ensembl{name}"


@attr.s()
class TranscriptInfo:
    so_rna_type = attr.ib(validator=is_a(str))
    regions: ty.List[SequenceRegion] = attr.ib(validator=is_a(list))
    from_gencode = attr.ib(validator=is_a(bool))

    @classmethod
    def from_feature(cls, record, feature) -> "TranscriptInfo":
        so_term = as_so_term(embl.rna_type(feature))
        return cls(
            so_rna_type=so_term,
            regions=helpers.regions(record, feature),
            from_gencode=False,
        )


@attr.s()
class FtpInfo:
    division = attr.ib(validator=is_a(Division))
    species = attr.ib(validator=is_a(str))
    data_files = attr.ib(validator=is_a(str))
    gff_file = attr.ib(validator=is_a(str))

    def writeable(self) -> ty.Tuple[str, str, str, str]:
        return (
            self.division.name,
            self.species,
            self.data_files,
            self.gff_file,
        )


@attr.s()
class Pseudogene:
    gene = attr.ib(validator=is_a(str))
    region = attr.ib(validator=is_a(SequenceRegion))

    def writeable(self) -> ty.Iterable[ty.List[str]]:
        return self.region.writeable(self.gene, is_upi=False)
