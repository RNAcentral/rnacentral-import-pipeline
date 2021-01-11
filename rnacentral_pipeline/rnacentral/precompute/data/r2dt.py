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

def maybe(fn, value):
    if value is None:
        return None
    return fn(value)


@attr.s(hash=True)
class R2dtHit:
    model_id = attr.ib(validator=is_a(int))
    model_name = attr.ib(validator=is_a(str))
    model_rna_type = attr.ib(validator=is_a(RnaType))
    sequence_coverage = attr.ib(validator=optional(is_a(float)))
    model_coverage = attr.ib(validator=optional(is_a(float)))
    sequence_basepairs = attr.ib(validator=optional(is_a(int)))
    model_basepairs = attr.ib(validator=optional(is_a(int)))

    @classmethod
    def build(cls, so_tree, raw):
        return cls(
            model_id=raw["model_id"],
            model_name=raw["model_name"],
            model_rna_type=RnaType.from_so_term(so_tree, raw["model_so_term"]),
            sequence_coverage=maybe(float, raw['sequence_coverage']),
            model_coverage=maybe(float, raw['model_coverage']),
            sequence_basepairs=maybe(int, raw['sequence_basepairs']),
            model_basepairs=maybe(int, raw['model_basepairs']),
        )

    def paired_ratio(self):
        if self.sequence_basepairs is None or self.model_basepairs is None:
            return None
        return float(self.sequence_basepairs) / float(self.model_basepairs)
