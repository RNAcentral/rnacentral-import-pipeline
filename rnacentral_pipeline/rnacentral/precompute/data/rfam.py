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
from attr.validators import optional

from rnacentral_pipeline.databases.data import RnaType


@attr.s(hash=True)
class HitComponent:
    completeness = attr.ib(validator=is_a(float), converter=float)
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))


@attr.s(hash=True)
class RfamHit:
    """
    This represents the information needed to represent an Rfam Hit to compute
    the QA information.
    """

    model = attr.ib(validator=is_a(str))
    model_rna_type = attr.ib(validator=is_a(RnaType))
    model_domain = attr.ib(validator=optional(is_a(str)))
    model_name = attr.ib(validator=is_a(str))
    model_long_name = attr.ib(validator=is_a(str))
    sequence_info = attr.ib(validator=is_a(HitComponent))
    model_info = attr.ib(validator=is_a(HitComponent))

    @classmethod
    def build(cls, so_tree, raw):
        """
        Create a new RfamHit object. This accepts a dict where all keys match
        the attributes of this class.
        """
        return cls(
            model=raw["model"],
            model_rna_type=RnaType.from_so_id(so_tree, raw["model_rna_type"]),
            model_domain=raw["model_domain"],
            model_name=raw["model_name"],
            model_long_name=raw["model_long_name"],
            sequence_info=HitComponent(
                completeness=raw["sequence_completeness"],
                start=raw["sequence_start"],
                stop=raw["sequence_stop"],
            ),
            model_info=HitComponent(
                completeness=raw["model_completeness"],
                start=raw["model_start"],
                stop=raw["model_stop"],
            ),
        )

    @property
    def url(self):
        return "http://rfam.org/family/%s" % self.model
