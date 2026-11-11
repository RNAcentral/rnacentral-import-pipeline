# -*- coding: utf-8 -*-

from __future__ import annotations

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
import typing as ty

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data import SequenceRegion
from rnacentral_pipeline.databases.data.utils import as_so_term
from rnacentral_pipeline.databases.ensembl import helpers
from rnacentral_pipeline.databases.helpers import embl

LOGGER = logging.getLogger(__name__)


@attr.s()
class TranscriptInfo:
    so_rna_type: str = attr.ib(validator=is_a(str))
    regions: ty.List[SequenceRegion] = attr.ib(validator=is_a(list))
    from_gencode: bool = attr.ib(validator=is_a(bool))

    @classmethod
    def from_feature(cls, record, feature) -> TranscriptInfo:
        rna_type = embl.rna_type(feature)
        if rna_type is None:
            raise ValueError(f"Did not get an RNA type for {feature}")
        so_term = as_so_term(rna_type)
        return cls(
            so_rna_type=so_term,
            regions=helpers.regions(record, feature),
            from_gencode=False,
        )


@attr.s()
class Pseudogene:
    gene = attr.ib(validator=is_a(str))
    region = attr.ib(validator=is_a(SequenceRegion))

    def writeable(self) -> ty.Iterable[ty.List[str]]:
        return self.region.writeable(self.gene, is_upi=False)
