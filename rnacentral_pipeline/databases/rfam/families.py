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

import re
import csv
import operator as op
import itertools as it

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .data import INFORMATIVE_NAMES
from .data import SO_TERM_MAPPING
from .data import RFAM_RNA_TYPE_MAPPING
from .data import DOMAIN_MAPPING


def empty_str_from(target):
    """
    Produces a function that will turn a string into an empty string if it
    matches the given target.
    """

    def func(raw):
        """
        The function which will do the conversion.
        """

        if raw == target:
            return None
        return raw
    return func


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    pretty_name = attr.ib(validator=is_a(str))
    so_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(str))
    domain = attr.ib()
    description = attr.ib(
        validator=optional(is_a(str)),
        converter=empty_str_from('NULL'),
    )
    seed_count = attr.ib(validator=is_a(int))
    full_count = attr.ib(validator=is_a(int))
    clan_id = attr.ib(converter=empty_str_from('NULL'))
    length = attr.ib(validator=is_a(int))

    @classmethod
    def from_dict(cls, row):
        family = row['id']
        domain = DOMAIN_MAPPING.get(family, None)

        clan = set(row['clan_id'].split(','))
        if len(clan) > 1:
            raise ValueError("Can only handle a single clan per family")
        clan = clan.pop()

        return cls(
            id=row['id'],
            name=row['name'],
            pretty_name=row['pretty_name'],
            so_terms=set(row['so_terms'].split(',')),
            rna_type=row['rna_type'].strip().strip(';'),
            domain=domain,
            description=row['description'],
            seed_count=int(row['seed_count']),
            full_count=int(row['full_count']),
            clan_id=clan,
            length=int(row['length']),
        )

    @property
    def is_suppressed(self):
        return 'lncRNA' in self.rna_type or \
            not self.rna_type.startswith('Gene')

    def guess_insdc_using_name(self):
        found = set()
        for name in INFORMATIVE_NAMES.keys():
            if re.search(name, self.name, re.IGNORECASE):
                found.add(name)

        if found:
            if len(found) > 1:
                raise ValueError("Name patterns not distinct %s" %
                                 ', '.join(sorted(found)))
            return INFORMATIVE_NAMES[found.pop()]
        return None

    def guess_insdc_using_rna_type(self):
        return RFAM_RNA_TYPE_MAPPING.get(self.rna_type, None)

    def guess_insdc_using_so_terms(self):
        terms = set(SO_TERM_MAPPING.get(so, None) for so in self.so_terms)
        if len(terms) > 1 and 'other' in terms:
            terms.remove('other')
        if len(terms) != 1:
            return None
        return terms.pop()

    def guess_insdc(self):
        return self.guess_insdc_using_name() or \
            self.guess_insdc_using_so_terms() or \
            self.guess_insdc_using_rna_type() or \
            'other'

    def writeable(self):
        return [
            self.id,
            self.name,
            self.pretty_name,
            self.description,
            self.clan_id,
            self.seed_count,
            self.full_count,
            self.length,
            self.domain,
            self.is_suppressed,
            self.guess_insdc(),
            self.rna_type,
        ]


def parse(handle):
    reader = csv.DictReader(handle, delimiter='\t')
    return map(RfamFamily.from_dict, reader)


def from_file(handle, output):
    data = parse(handle)
    data = map(op.methodcaller('writeable'), data)
    writer = csv.writer(output)
    writer.writerows(data)
