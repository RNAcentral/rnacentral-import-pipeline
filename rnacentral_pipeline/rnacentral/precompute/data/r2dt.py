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

from rnacentral_pipeline.databases.data import RnaType


@attr.s()
class R2dtHit:
    model_id = attr.ib(validator=is_a(int))
    model_name = attr.ib(validator=is_a(str))
    model_rna_type = attr.ib(validator=is_a(RnaType))

    def build(cls, raw):
        return cls(
            model_id=raw["model_id"],
            model_name=raw["model_name"],
            model_rna_type=RnaType.from_so_term(raw["model_so_term"]),
        )
