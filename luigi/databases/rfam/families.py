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
import itertools as it

import attr
from attr.validators import optional
from attr.validators import instance_of as is_a

from .data import INFORMATIVE_NAMES
from .data import SO_TERM_MAPPING
from .data import RFAM_RNA_TYPE_MAPPING
from .data import DOMAIN_MAPPING

QUERY = """
select
    family.rfam_acc id,
    family.rfam_id as name,
    family.description as pretty_name,
    group_concat(concat('SO:', dbs.db_link)) as so_terms,
    family.type as rna_type,
    family.comment as description,
    family.num_seed as seed_count,
    family.num_full as full_count,
    group_concat(mem.clan_acc) as clan_id,
    family.clen as length
from family
left join database_link dbs on dbs.rfam_acc = family.rfam_acc
left join clan_membership mem on mem.rfam_acc = family.rfam_acc
where
    dbs.db_id = 'SO'
group by family.rfam_acc
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
            return None
        return raw
    return func


@attr.s(frozen=True)
class RfamFamily(object):
    id = attr.ib(validator=is_a(basestring))
    name = attr.ib(validator=is_a(basestring))
    pretty_name = attr.ib(validator=is_a(basestring))
    so_terms = attr.ib(validator=is_a(set))
    rna_type = attr.ib(validator=is_a(basestring))
    domain = attr.ib()
    description = attr.ib(
        validator=optional(is_a(basestring)),
        convert=empty_str_from('NULL'),
    )
    seed_count = attr.ib(validator=is_a(int))
    full_count = attr.ib(validator=is_a(int))
    clan_id = attr.ib(convert=empty_str_from('NULL'))
    length = attr.ib(validator=is_a(int))

    @classmethod
    def from_dict(cls, row):
        family = row['id']
        domain = DOMAIN_MAPPING.get(family, None)
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
            clan_id=row['clan_id'],
            length=int(row['length']),
        )

    @property
    def is_suppressed(self):
        return 'lncRNA' in self.rna_type or \
            not self.rna_type.startswith('Gene')

    def guess_insdc_using_name(self):
        found = set()
        for name in INFORMATIVE_NAMES.iterkeys():
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

    def writeable(self):
        data = attr.asdict(self)
        del data['so_terms']
        data['short_name'] = data.pop('name')
        data['long_name'] = data.pop('pretty_name')
        data['is_suppressed'] = int(self.is_suppressed)
        data['rfam_rna_type'] = data['rna_type']
        data['rna_type'] = self.guess_insdc()
        yield data


def parse(handle):
    reader = csv.DictReader(handle, delimiter='\t')
    return it.imap(RfamFamily.from_dict, reader)
