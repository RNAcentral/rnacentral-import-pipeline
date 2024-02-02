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

import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


@attr.s(frozen=True)
class OntologyTerm:
    """
    This represents a single term in a specific ontology.
    """

    ontology: str = attr.ib(validator=is_a(str), converter=str)
    ontology_id: str = attr.ib(validator=is_a(str), converter=str)
    name: str = attr.ib(validator=is_a(str), converter=str)
    definition: ty.Optional[str] = attr.ib(validator=optional(is_a(str)))
    synonyms: ty.List[str] = attr.ib(validator=is_a(list))
    insdc_qualifier: ty.Optional[str] = attr.ib(
        validator=optional(is_a(str)),
        default=None,
    )

    def writeable(self):
        return [
            self.ontology_id,
            self.ontology,
            self.name,
            self.definition,
        ]

    def writeables(self):
        yield self.writeable()
