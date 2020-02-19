# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import re
import enum
import typing as ty
from datetime import date

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .references import IdReference


@enum.unique
class InteractorType(enum.Enum):
    A = enum.auto()
    B = enum.auto()


@attr.s(frozen=True)
class InteractionIdentifier:
    key = attr.ib(validator=is_a(str))
    value = attr.ib(validator=is_a(str))
    name = attr.ib(validator=optional(is_a(str)))

    def full_id(self) -> str:
        name = ''
        if self.name:
            name = ' (%s)' % self.name
        return "%s%s" % (self.unique_id(), name)

    def unique_id(self) -> str:
        quote = ''
        if ':' in self.value:
            quote = '"'

        return '{key}:{quote}{value}{quote}'.format(
            key=self.key,
            value=self.value,
            quote=quote,
        )


@attr.s(frozen=True)
class Interactor:
    id = attr.ib(validator=is_a(InteractionIdentifier))
    alt_ids: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    aliases: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    taxid = attr.ib(validator=optional(is_a(int)))
    biological_role: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    experimental_role: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    interactor_type: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    xrefs: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    annotations: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    features: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    stoichiometry = attr.ib(validator=optional(is_a(int)))
    participant_identification: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))

    @property
    def database(self):
        return self.id.key

    @property
    def external_id(self):
        return self.id.value


@attr.s(frozen=True)
class Interaction:
    ids: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    interactor1 = attr.ib(validator=is_a(Interactor))
    interactor2 = attr.ib(validator=is_a(Interactor))
    methods: ty.List[str] = attr.ib(validator=is_a(list))
    types: ty.List[str] = attr.ib(validator=is_a(list))
    xrefs: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    annotations: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    confidence: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    source_database: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    is_negative = attr.ib(validator=is_a(bool))
    publications: ty.List[IdReference] = attr.ib(validator=is_a(list))
    create_date = attr.ib(validator=is_a(date))
    update_date = attr.ib(validator=is_a(date))
    host_organisms = attr.ib(validator=is_a(int))

    def writeables(self):
        yield []
