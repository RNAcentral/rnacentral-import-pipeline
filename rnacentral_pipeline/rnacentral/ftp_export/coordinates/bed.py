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

import six

from rnacentral_pipeline.databases.data import regions

from . import data as coord


@attr.s(slots=True, frozen=True)
class BedEntry(object):
    rna_id = attr.ib(validator=is_a(str))
    rna_type = attr.ib(validator=is_a(str))
    databases = attr.ib(validator=is_a(str))
    region = attr.ib(validator=is_a(regions.SequenceRegion))
    score = attr.ib(default=0, validator=is_a(int))
    rgb = attr.ib(default=(63, 125, 151), validator=is_a(tuple))

    @classmethod
    def from_coordinate(cls, coordinate):
        return cls(
            rna_id=coordinate.rna_id,
            rna_type=coordinate.metadata['rna_type'],
            databases=','.join(coordinate.metadata['databases']),
            region=coordinate.region.as_zero_based(),
        )

    @property
    def bed_chromosome(self):
        if self.region.chromosome in ['MT', 'chrMT']:
            return 'chrM'
        return 'chr' + self.region.chromosome

    @property
    def bed_rgb(self):
        return ','.join(str(c) for c in self.rgb)

    def sizes(self):
        return self.region.sizes()

    def starts(self):
        starts = []
        end = self.region.exons[0].start
        for exon in self.region.exons[1:]:
            start = exon.start - end
            assert start > 0, "Invalid start for %s" % (self)
            starts.append(start)
        return [0] + starts

    def writeable(self):
        return [
            self.bed_chromosome,
            self.region.start,
            self.region.stop,
            self.rna_id,
            self.score,
            self.region.strand.display_string(),
            self.region.start,
            self.region.stop,
            self.bed_rgb,
            len(self.region.exons),
            ','.join(str(s) for s in self.sizes()),
            ','.join(str(s) for s in self.starts()),
            '.',
            self.rna_type,
            self.databases,
        ]


def from_json(handle, out):
    """
    Transform raw coordinate data into bed format.
    """

    data = coord.from_file(handle)
    data = six.moves.map(BedEntry.from_coordinate, data)
    data = six.moves.map(op.methodcaller('writeable'), data)
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerows(data)
