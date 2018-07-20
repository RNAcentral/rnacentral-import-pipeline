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

from . import data as coord


@attr.s(frozen=True)
class BedBlock(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def from_endpoint(cls, endpoint):
        return cls(
            start=endpoint.start,
            stop=endpoint.stop,
        )

    @property
    def size(self):
        start = self.start - 1
        return (self.stop - start) or 1


@attr.s(slots=True, frozen=True)
class BedEntry(object):
    chromosome = attr.ib(validator=is_a(basestring))
    rna_id = attr.ib(validator=is_a(basestring))
    blocks = attr.ib(validator=is_a(list))
    strand = attr.ib(validator=is_a(int))
    rna_type = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    coord_system = attr.ib(validator=is_a(basestring))
    score = attr.ib(default=0, validator=is_a(int))

    @classmethod
    def from_region(cls, region):
        return cls(
            chromosome=region.chromosome,
            rna_id=region.rna_id,
            blocks=[BedBlock.from_endpoint(e) for e in region.endpoints],
            strand=region.strand,
            rna_type=region.rna_type,
            databases=region.databases,
            coord_system=region.coord_system,
        )

    def block_sizes(self):
        return [b.size for b in self.blocks]

    def block_starts(self):
        return [0] + [b.start - self.start for b in self.blocks[1:]]

    def bed_block_sizes(self):
        return ','.join(self.block_sizes())

    def bed_block_starts(self):
        return ','.join(self.block_starts())

    @property
    def start(self):
        return self.blocks[0].start

    @property
    def stop(self):
        return self.blocks[-1].stop

    def bed_strand(self):
        if self.strand == 1:
            return '+'
        if self.strand == -1:
            return '-'
        raise ValueError('Unknown Strand')

    @property
    def bed_chromosome(self):
        chromosome = self.chromosome
        if self.coord_system == 'chromosome':
            chromosome = 'chr' + chromosome
        if chromosome in ['MT', 'chrMT']:
            chromosome = 'chrM'
        return chromosome

    def writeable(self):
        return [
            self.bed_chromosome,
            self.start,
            self.stop,
            self.rna_id,
            self.score,
            self.bed_strand,
            self.start,
            self.stop,
            '63,125,151',
            len(self.blocks),
            self.bed_block_sizes(),
            self.bed_block_starts(),
            '.',
            self.rna_type,
            self.databases,
        ]


def located_sequences_as_bed(sequences):
    """
    Transform the iterable of LocatedSequence into an iterable of BedEntry
    objects.
    """

    for sequence in sequences:
        for region in sequence.regions:
            yield BedEntry.from_region(region)


def from_json(handle, out):
    """
    Transform raw coordinate data into bed format.
    """

    data = coord.from_file(handle)
    bed = located_sequences_as_bed(data)
    bed = it.imap(op.methodcaller('writeable'), bed)
    writer = csv.writer(out, delimiter='\t')
    writer.writerows(bed)
