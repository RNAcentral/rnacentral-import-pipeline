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
import json
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

    def from_rnacentral(self) -> bool:
        return self.key == 'rnacentral'

    def simple_id(self) -> str:
        return '%s:%s' % (self.key, self.value)


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
    annotations = attr.ib(validator=is_a(str))
    features: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))
    stoichiometry = attr.ib(validator=optional(is_a(int)))
    participant_identification: ty.List[InteractionIdentifier] = attr.ib(validator=is_a(list))

    @property
    def database(self):
        return self.id.key

    @property
    def external_id(self):
        return self.id.value

    def from_rnacentral(self) -> bool:
        if self.id.from_rnacentral():
            return True
        return any(a.from_rnacentral() for a in self.alt_ids) or \
            any(a.from_rnacentral() for a in self.aliases)

    def all_ids(self):
        return [self.id] + self.alt_ids + self.aliases

    def urs_taxid(self) -> str:
        found = set()
        for id in self.all_ids():
            if id.from_rnacentral():
                found.add(id.value)

        if len(found) == 1:
            return found.pop()

        if len(found) > 1:
            raise ValueError("Found more than one URS/taxid")

        raise ValueError("Could not find a URS/taxid")

    def names(self) -> ty.List[str]:
        names: ty.Set[str] = set()
        names.update(id.value for id in self.all_ids())
        return sorted(names)


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

    def involves_rnacentral(self) -> bool:
        return self.interactor1.from_rnacentral() or \
            self.interactor2.from_rnacentral()

    def intact_id(self):
        for id in self.ids:
            if id.key == 'intact':
                return id
        return None

    @property
    def rnacentral_iteractor(self):
        if self.interactor1.from_rnacentral():
            return self.interactor1
        elif self.interactor2.from_rnacentral():
            return self.interactor2
        raise ValueError("No rnacentral iteractor")

    @property
    def urs_taxid(self):
        return self.rnacentral_iteractor.urs_taxid()

    def writeable(self):
        other = None
        urs_taxid = None
        if self.interactor1.from_rnacentral():
            urs_taxid = self.interactor1.urs_taxid()
            other = self.interactor2
        elif self.interactor2.from_rnacentral():
            urs_taxid = self.interactor2.urs_taxid()
            other = self.interactor1
        else:
            raise ValueError("Cannot write iteraction")

        intact_id = self.intact_id()
        if intact_id:
            intact_id = intact_id.value
        return [
            intact_id,
            urs_taxid,
            other.id.simple_id(),
            json.dumps(other.names()),
            other.taxid,
        ]
