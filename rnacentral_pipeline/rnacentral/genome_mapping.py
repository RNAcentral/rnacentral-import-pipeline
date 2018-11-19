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


import csv
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a

from Bio import SearchIO

from rnacentral_pipeline.databases.data import regions


@attr.s()
class Hit(object):
    chromsome = attr.ib(validator=is_a(basestring))
    upi = attr.ib(validator=is_a(basestring))
    sequence_length = attr.ib(validator=is_a(int))
    matches = attr.ib(validator=is_a(int))
    target_insertions = attr.ib(validator=is_a(int))
    strand = attr.ib(validator=is_a(int))

    @classmethod
    def from_biopython(cls, hit):
        pass

    @property
    def region_id(self):
        return regions.region_id(self)

    @property
    def exons(self):
        return []

    def writeable(self):
        return regions.write_locations(self, None)


def select_possible(hit):
    if hit.matches < 100 and hit.target_insertions > 10:
        return False
    if hit.matches == hit.sequence_length:
        return True
    if hit.sequence_length > 15 and hit.match_fraction > 0.95:
        return True
    return False


def select_best(hits):
    best = max(hits, key=op.attrgetter('match_fraction'))
    if best.match_fraction == 1.0:
        return [h for h in hits if h.match_fraction == best]
    return hits


def select_hits(handle):
    query_results = SearchIO.parse(handle, 'blat-psl')
    hits = it.imap(Hit.from_biopython, query_results)
    hits = it.ifilter(select_possible, hits)
    hits = it.groupby(hits, op.attrgetter('upi'))
    hits = it.imap(op.itemgetter(1), hits)
    hits = it.imap(select_best, hits)
    hits = it.chain.from_iterable(hits)
    return hits


def write_selected(hits, output):
    selected = select_hits(hits)
    selected = it.imap(op.methodcaller('writeable'), selected)
    selected = it.chain.from_iterable(selected)
    csv.writer(output).writerows(selected)
