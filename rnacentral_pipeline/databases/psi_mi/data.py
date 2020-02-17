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
import typing as ty

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

IDENTIFIER_PATTERN = re.compile(r'^(.+?):"?(.+?)"?(\((.+?)\))?$')


@attr.s()
class InteractionAnnotation:
    pass


@attr.s()
class Identifier:
    key = attr.ib(validator=is_a(str))
    value = attr.ib(validator=is_a(str))
    name = attr.ib(validator=optional(is_a(str)))

    @classmethod
    def build(cls, raw):
        matches = re.match(IDENTIFIER_PATTERN, raw)
        if not matches:
            raise ValueError("Invalid identifier format: %s" % raw)

        name = matches.group(3)
        if name:
            name = name.replace('(', '').replace(')', '')

        return cls(
            key=matches.group(1), 
            value=matches.group(2),
            name=name,
        )

    def full_id(self):
        name = ''
        if self.name:
            name = ' (%s)' % self.name
        return "%s%s" % (self.unique_id(), name)

    def unique_id(self):
        quote = ''
        if ':' in self.value:
            quote = '"'
        return '{key}:{quote}{value}{quote}'.format(
            key=self.key,
            value=self.value,
            quote=quote,
        )


@attr.s()
class Interactor:
    id: str = attr.ib(validator=is_a(str))
    alt_ids: ty.List[str] = attr.ib(validator=is_a(list))
    aliases: ty.List[str] = attr.ib(validator=is_a(list))
    taxid: int = attr.ib(validator=is_a(int))
    biological_role: str = attr.ib()
    experimental_role: str = attr.ib()
    interactor_type: str = attr.ib()
    xrefs: str = attr.ib()
    annotations: str = attr.ib()
    checksum: str = attr.ib()
    features: str = attr.ib()
    stoichiometry: str = attr.ib()
    participant_identification: str = attr.ib()
    biological_effect: str = attr.ib()

    @classmethod
    def from_tab(cls, row, offset):
        parts = {
            'unique_ids': row[offset],
            'alternatives': row[offset + 2],
            'aliases': row[offset + 4],
            'taxids': row[offset + 9],
            'xrefs': row[offset + 22]
        }

        unique_ids = identifiers(parts['unique_ids'])
        assert len(unique_ids) == 1

        taxids = identifiers(parts['taxid'])
        assert len(taxids) == 1

        return cls(
            id=unique_ids[0],
            taxid=taxids[0].value,
            alt_ids=identifiers(parts['identifiers']),
            aliases=identifiers(parts['aliases']),
            xrefs=identifiers(parts['xrefs']),
        )

    @classmethod
    def database(self):
        return self.id.key

    @property
    def external_id(self):
        return self.id.value


@attr.s()
class SourceDatabase:
    id: str = attr.ib(validator=is_a(str))
    name: str = attr.ib(validator=is_a(str))

    @classmethod
    def build(cls, raw):
        return cls(


@attr.s()
class Interaction:
    ids: ty.List[str] = attr.ib()
    interactor1 = attr.ib(validator=is_a(Interactor))
    interactor2 = attr.ib(validator=is_a(Interactor))
    methods: ty.List[str] = attr.ib(validator=is_a(list))
    reference: ty.List[str] = attr.ib(validator=is_a(list))
    interaction_types: ty.List[str] = attr.ib()
    source_database: ty.List[SourceDatabase] = attr.ib()
    is_negative = attr.ib(validator=is_a(bool))
    publications = attr.ib(validator=is_a(list))

    @classmethod
    def from_tab(cls, row):
        return cls(
            ids=[],
            interactor1=Interactor.from_tab(row, 1),
            interactor1=Interactor.from_tab(row, 2),
            methods=[],
            reference=[],
            source_database=SourceDatabase.build(row
        )
