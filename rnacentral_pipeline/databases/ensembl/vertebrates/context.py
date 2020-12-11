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

import gffutils
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.rfam import utils as rfutils

from rnacentral_pipeline.databases.ensembl.gff import load_coordinates
from rnacentral_pipeline.databases.ensembl.vertebrates.rna_type_inference import RnaTypeInference

LOGGER = logging.getLogger(__name__)


@attr.s(frozen=True, slots=True)
class Context:
    supressed_mapping = attr.ib()
    inference = attr.ib(validator=is_a(RnaTypeInference))
    rfam_names: ty.Dict[str, str] = attr.ib(validator=is_a(dict))
    excluded = attr.ib(validator=is_a(set))
    gff = attr.ib(validator=is_a(SqliteDict))

    @classmethod
    def build(cls, gff_file, family_file=None, excluded_file=None):
        """
        Create a Context, by parsing the family file (families from Rfam), and
        the gencode_file (gff3 file from GENCODE) and excluded_file a list of
        Ensembl Ids to ignore.
        """

        if not family_file:
            raise ValueError("Must give family file for E!vert")

        with open(family_file, 'r') as raw:
            supressed_mapping = rfutils.name_to_suppression(raw)

        with open(family_file, 'r') as raw:
            supressed_mapping.update(rfutils.id_to_suppression(raw))

        with open(family_file, 'r') as raw:
            inference = RnaTypeInference(raw)

        with open(family_file, 'r') as raw:
            rfam_names = rfutils.id_to_pretty_name(raw)

        with open(family_file, 'r') as raw:
            rfam_names.update(rfutils.name_to_pretty_name(raw))

        excluded = set()
        if excluded_file:
            excluded = set(l.strip() for l in excluded_file)

        return cls(
            supressed_mapping,
            inference,
            rfam_names=rfam_names,
            excluded=excluded,
            gff=load_coordinates(gff_file),
        )

    def rfam_name(self, locus_tag: str, default=None):
        return self.rfam_names.get(locus_tag, default)

    def is_supressed(self, feature) -> bool:
        rfam_models = self.inference.rfam_xref(feature)
        if not rfam_models:
            return False

        for rfam_model in rfam_models:
            name = self.inference.rfam_name(rfam_model)
            if name is None:
                continue
            if name not in self.supressed_mapping:
                raise ValueError("Unknown Rfam model name: %s", name)
            if self.supressed_mapping[name]:
                return True
        return False

    def from_gencode(self, entry: data.Entry) -> bool:
        return self.gff[entry.primary_id].from_gencode

    def is_excluded(self, entry: data.Entry) -> bool:
        return entry.accession in self.excluded
