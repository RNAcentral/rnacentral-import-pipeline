# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

import six

from rnacentral_pipeline.databases.data import IdReference

IDENTIFIER_PATTERN = re.compile(r'^(.+?):"?(.+?)"?(?\((.+?)\))?$')


class Identifier(object):
    key = attr.ib(validator=is_a(six.text_type))
    value = attr.ib(validator=is_a(six.text_type))
    name  = attr.ib(validator=optional(is_a(six.text_type)))

    def full_id(self):
        pass

    def unique_id(self):
        quote = ''
        if ':' in self.value:
            quote = '"'
        return '{key}:{quote}{value}{quote}'.format(
            key=self.key,
            value=self.value,
            quote=quote,
        )

    @classmethod
    def bulid(cls, raw):
        matches = re.match(IDENTIFIER_PATTERN, raw)
        if not matches:
            raise ValueError("Invalid identifier format: %s" % raw)
        return cls(
            key=matches.group(1), 
            value=matches.group(2),
            name=matches.group(3),
        )


class Interactor(object):
    identifier = attr.ib(validator=is_a(Identifier))
    taxonomy_id = attr.ib(validator=is_a(six.integer_types))
    alternatives = attr.ib(validator=is_a(list))
    aliases = attr.ib(validator=is_a(list))
    xrefs = attr.ib(validator=is_a(list))

    @property
    def database(self):
        return self.identifier.key

    @property
    def external_id(self):
        return self.identifier.value

    def unique_id(self):
        return '%s:%s' % (self.database, self.external_id)


@attr.s()
class Interaction(object):
    interactor = attr.ib(validator=is_a(Interactor))
    urs = attr.ib(validator=is_a(Interactor))
    interaction = attr.ib(validator=is_a(six.text_type))
    publications = attr.ib(validator=is_a(list))
    is_negative = attr.ib(validator=is_a(bool))

    def writeables(self):
        for publication in self.publications:
            yield [
                self.urs.unique_id(),
                self.interactor.unique_id(),
                self.interaction,
                publication.normalized_id,
                self.is_negative,
            ]

    def writeable_ref_ids(self):
        for publication in self.publications:
            yield publication.normalized_id

    def writeable_term(self):
        return [interaction.external_id]
