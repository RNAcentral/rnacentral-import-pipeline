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
import collections as coll

import attr
from attr.validators import instance_of as is_a

from .data import INFORMATIVE_NAMES
from .data import SO_TERM_MAPPING
from .data import RFAM_RNA_TYPE_MAPPING
from .data import DOMAIN_MAPPING

QUERY = """
"""


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
            return u''
        return raw
    return func


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(basestring))
    name = attr.ib(validator=is_a(basestring))
    pretty_name = attr.ib(validator=is_a(basestring))
    so_terms = attr.ib(validator=is_a(set))
    go_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(basestring))
    domain = attr.ib()
    description = attr.ib(
        validator=is_a(basestring),
        convert=empty_str_from(r'\N'),
    )
    seed_count = attr.ib(validator=is_a(int))
    full_count = attr.ib(validator=is_a(int))
    clan_id = attr.ib()
    length = attr.ib(validator=is_a(int))

    @classmethod
    def from_row(cls, row):
        pass

    @classmethod
    def build_all(cls, clan_file, link_file, family_file):
        so_terms = coll.defaultdict(set)
        go_terms = coll.defaultdict(set)
        for line in link_file:
            parts = line.split()
            if parts[1] == 'SO':
                so_terms[parts[0]].add('SO:%s' % parts[2])
            if parts[1] == 'GO':
                go_name = ' '.join(parts[3:])
                go_terms[parts[0]].add(('GO:%s' % parts[2], go_name))

        clans = coll.defaultdict(set)
        for line in clan_file:
            parts = line.split()
            clans[parts[1]] = parts[0]

        families = []
        for row in csv.reader(family_file, delimiter='\t'):
            family = row[0]
            domain = DOMAIN_MAPPING.get(family, None)
            families.append(cls(
                id=family,
                name=row[1],
                pretty_name=row[3],
                so_terms=so_terms[family],
                go_terms=go_terms[family],
                rna_type=row[18].strip().strip(';'),
                domain=domain,
                description=row[9],
                seed_count=int(row[14]),
                full_count=int(row[15]),
                clan_id=clans.get(family, None),
                length=int(row[28]),
            ))
        return families

    @property
    def is_suppressed(self):
        return 'lncRNA' in self.rna_type or \
            not self.rna_type.startswith('Gene')

    def guess_insdc_using_name(self):
        found = set()
        for name, _ in INFORMATIVE_NAMES.items():
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
            self.guess_insdc_using_rna_type()


def build_all(cls, clan_file, membership_file):
    membership = coll.defaultdict(set)
    for line in membership_file:
        parts = line.split()
        membership[parts[0]].add(parts[1])
    membership = dict(membership)

    clans = []
    for line in clan_file:
        parts = line.strip().split('\t')
        clans.append(cls(
            id=parts[0],
            name=parts[3],
            description=parts[5],
            families=membership[parts[0]],
        ))
    return clans


def parse(handle):
    reader = csv.DictReader(handle)
    for row in reader:
        RfamFamily.from_dict(row)

    # family_file = get_family_file(version=version)
    # link_file = get_link_file(version=version)
    # clan_file = get_clan_membership(version=version)
    return RfamFamily.build_all(clan_file, link_file, family_file)
