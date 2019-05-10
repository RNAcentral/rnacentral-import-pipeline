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
import logging

import six

import attr
from attr.validators import instance_of as is_a

from rnacentral_pipeline.databases.data.regions import Exon
from rnacentral_pipeline.databases.data.regions import Strand
from rnacentral_pipeline.databases.data.regions import SequenceRegion
from rnacentral_pipeline.databases.data.regions import CoordinateSystem

LOGGER = logging.getLogger(__name__)

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
class BlatHit(object):
    upi = attr.ib(validator=is_a(six.text_type), converter=six.text_type)
    sequence_length = attr.ib(validator=is_a(int))
    matches = attr.ib(validator=is_a(int))
    target_insertions = attr.ib(validator=is_a(int))
    region = attr.ib(validator=is_a(SequenceRegion))

    @classmethod
    def build(cls, assembly_id, raw):
        parts = six.moves.zip(raw['tStarts'], raw['blockSizes'])
        exons = [Exon(s, s+l) for (s, l) in parts]

        return cls(
            upi=raw['qName'],
            sequence_length=raw['qSize'],
            matches=raw['matches'],
            target_insertions=raw['tBaseInsert'],
            region=SequenceRegion(
                assembly_id=assembly_id,
                chromosome=raw['tName'],
                strand=raw['strand'],
                exons=exons,
                coordinate_system=CoordinateSystem.zero_based(),
            ),
        )

    @property
    def name(self):
        return self.region.name(self.upi, is_upi=True)

    @property
    def match_fraction(self):
        return float(self.matches) / float(self.sequence_length)

    def writeable(self):
        return self.region.writeable(self.upi, is_upi=True)


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


def parse_psl(assembly_id, handle):
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

        yield BlatHit.build(assembly_id, result)


def parse_json(handle):
    for line in handle:
        data = json.loads(line)
        yield BlatHit(**data)


def select_hits(hits):
    hits = six.moves.filter(select_possible, hits)
    hits = it.groupby(hits, op.attrgetter('upi'))
    hits = six.moves.map(op.itemgetter(1), hits)
    hits = six.moves.map(select_best, hits)
    hits = it.chain.from_iterable(hits)
    return hits


def write_json(hits, output):
    for hit in hits:
        output.write(json.dumps(hit))
        output.write('\n')


def write_importable(hits, output):
    writeable = six.moves.map(op.methodcaller('writeable'), hits)
    writeable = it.chain.from_iterable(hits)
    csv.writer(output).writerows(writeable)


def as_json(assembly_id, hits, output):
    parsed = parse_psl(assembly_id, hits)
    parsed = six.moves.map(attr.asdict, selected)
    write_json(parsed, output)


def select_json(hits, output, sort=False):
    parsed = parse_json(hits)
    if sort:
        parsed = sorted(parsed, key=op.itemgetter('upi'))
    selected = select_hits(parsed)
    write_json(selected, output)
