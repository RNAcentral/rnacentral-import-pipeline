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

import logging

import attr
from attr.validators import instance_of as is_a

import gffutils

from rnacentral_pipeline.databases.rfam import utils as rfutils

from .rna_type_inference import RnaTypeInference

LOGGER = logging.getLogger(__name__)


@attr.s(frozen=True, slots=True)
class Context(object):
    supressed_mapping = attr.ib()
    inference = attr.ib()
    gencode_ids = attr.ib(validator=is_a(set))
    rfam_names = attr.ib(validator=is_a(dict))
    excluded = attr.ib(validator=is_a(set))

    @classmethod
    def build(cls, family_file, gencode_file=None, excluded_file=None):
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

        gencode_ids = set()
        if gencode_file:
            gff = gffutils.create_db(gencode_file, ':memory:',
                                     merge_strategy='warning')
            gencode_ids = {f['ID'][0] for f in gff.features_of_type('transcript')}

        excluded = set()
        if excluded_file:
            with open(excluded_file, 'r') as raw:
                excluded = set(l.strip() for l in raw)

        return cls(
            supressed_mapping,
            inference,
            gencode_ids=gencode_ids,
            rfam_names=rfam_names,
            excluded=excluded,
        )

    def rfam_name(self, locus_tag, default=None):
        return self.rfam_names.get(locus_tag, default)

    def is_supressed(self, feature):
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

    def from_gencode(self, entry):
        return entry.accession in self.gencode_ids

    def is_excluded(self, entry):
        return entry.accession in self.excluded
