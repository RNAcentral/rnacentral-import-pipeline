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
import tempfile
import logging
import typing as ty

import attr
from attr.validators import instance_of as is_a

from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.rfam import utils as rfutils

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
        return f'Ensembl{name}'


@attr.s()
class TranscriptInfo:
    so_rna_type = attr.ib(validator=is_a(str))
    regions: ty.List[data.SequenceRegion] = attr.ib(validator=is_a(list))
    from_gencode = attr.ib(validator=is_a(bool))
