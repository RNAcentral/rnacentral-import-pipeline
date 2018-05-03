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

import csv
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a

QUERY = """
select
    clan.*,
    membership.rfam_acc
from clan
join clan_membership membership
on
    membership.clan_acc = clan.clan_acc
order by clan.clan_acc
"""


@attr.s()
class RfamClan(object):
    """
    A class to represent the information about Rfam clans.
    """

    id = attr.ib(validator=is_a(str))
    name = attr.ib(validator=is_a(str))
    description = attr.ib(validator=is_a(str))
    families = attr.ib(validator=is_a(set))

    @property
    def family_count(self):
        """
        The number of families in this clan.
        """

        return len(self.families)

    def writeable(self):
        data = attr.asdict(self)
        data['family_count'] = clan.family_count
        return data


def parse(handle):
    """
    Parse the given file and produce an iterable of all Rfam clans in that
    file. The file should be a tsv file as a result of running QUERY.
    """

    reader = csv.DictReader(handle, delimiter='\t')
    grouped = it.groupby(op.itemgetter('clan_acc'), reader)
    for clan_acc, entries in grouped:
        entries = list(entries)
        yield RfamClan(
            id=clan_acc,
            name=entries[0]['name'],
            description=entries[0]['description'],
            families=[e['rfam_acc'] for e in entries],
        )
