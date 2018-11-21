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

from rnacentral_pipeline.databases.data import regions

FIELDS = [
    'matches',  # Number of bases that match that aren't repeats
    'misMatches',  # Number of bases that don't match
    'repMatches',  # Number of bases that match but are part of repeats
    'nCount',  # Number of "N" bases
    'qNumInsert',  # Number of inserts in query
    'qBaseInsert',  # Number of bases inserted in query
    'tNumInsert',  # Number of inserts in target
    'tBaseInsert',  # Number of bases inserted in target
    'strand',  # "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
    'qName',  # Query sequence name
    'qSize',  # Query sequence size.
    'qStart',  # Alignment start position in query
    'qEnd',  # Alignment end position in query
    'tName',  # Target sequence name
    'tSize',  # Target sequence size
    'tStart',  # Alignment start position in target
    'tEnd',  # Alignment end position in target
    'blockCount',  # Number of blocks in the alignment (a block contains no gaps)
    'blockSizes',  # Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
    'qStarts',  # Comma-separated list of starting positions of each block in query
    'tStarts',  # Comma-separated list of starting positions of each block in target
]


@attr.s()
class Hit(object):
    assembly_id = attr.ib(validator=is_a(basestring))
    chromosome = attr.ib(validator=is_a(basestring))
    upi = attr.ib(validator=is_a(basestring))
    sequence_length = attr.ib(validator=is_a(int))
    matches = attr.ib(validator=is_a(int))
    target_insertions = attr.ib(validator=is_a(int))
    strand = attr.ib(validator=is_a(int))
    exons = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, assembly_id, raw):
        exons = []
        for (start, size) in zip(raw['tStarts'], raw['blockSizes']):
            exons.append(regions.Exon(
                start=start + 1,
                stop=start+size,
            ))

        return cls(
            assembly_id=assembly_id,
            chromosome=raw['tName'],
            upi=raw['qName'],
            sequence_length=raw['qSize'],
            matches=raw['matches'],
            target_insertions=raw['tBaseInsert'],
            strand=raw['strand'],
            exons=exons,
        )

    @property
    def name(self):
        return regions.region_id(self)

    @property
    def match_fraction(self):
        return float(self.matches) / float(self.sequence_length)

    def writeable(self):
        return regions.write_locations(self, None)


def select_possible(hit):
    if hit.matches < 100 and hit.target_insertions > 25:
        return False
    if hit.matches == hit.sequence_length:
        return True
    if hit.sequence_length > 15 and hit.match_fraction > 0.95 and \
            hit.match_fraction < 1:
        return True
    return False


def select_best(hits):
    hits = list(hits)
    best = max(hits, key=op.attrgetter('match_fraction'))
    if best.match_fraction == 1.0:
        return [h for h in hits if h.match_fraction == best.match_fraction]
    return hits


def parse(assembly_id, handle):
    to_split = ['blockSizes', 'qStarts', 'tStarts']
    for row in csv.reader(handle, delimiter='\t'):
        result = dict(zip(FIELDS, row))
        for key in to_split:
            result[key] = [int(v) for v in result[key].split(',') if v]
        lens = {len(result[k]) for k in to_split}
        assert len(lens) == 1

        for key, value in result.items():
            if key not in to_split and 'Name' not in key and key != 'strand':
                result[key] = int(value)

        result['strand'] = regions.as_strand(result['strand'])
        yield Hit.build(assembly_id, result)


def select_hits(assembly_id, handle):
    hits = parse(assembly_id, handle)
    hits = it.ifilter(select_possible, hits)
    hits = it.groupby(hits, op.attrgetter('upi'))
    hits = it.imap(op.itemgetter(1), hits)
    hits = it.imap(select_best, hits)
    hits = it.chain.from_iterable(hits)
    return hits


def write_selected(assembly_id, hits, output):
    selected = select_hits(assembly_id, hits)
    selected = it.imap(op.methodcaller('writeable'), selected)
    selected = it.chain.from_iterable(selected)
    csv.writer(output).writerows(selected)
