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

import attr
from attr.validators import instance_of as is_a


class UnknownStrand(Exception):
    """
    Raised when a strand integer has an invalid value.
    """
    pass


def as_strand(value):
    if isinstance(value, int):
        return value
    elif isinstance(value, float):
        return int(value)
    elif isinstance(value, basestring):
        if value == '+' or value == '1':
            return 1
        elif value == '-' or value == '-1':
            return -1
    raise UnknownStrand("No way to handle strand: " + str(value))


@attr.s(frozen=True)
class Exon(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))


@attr.s()
class SequenceRegion(object):
    chromosome = attr.ib(validator=is_a(basestring))
    strand = attr.ib(validator=is_a(int), convert=as_strand)
    exons = attr.ib(validator=is_a(list))
    assembly_id = attr.ib(validator=is_a(basestring))

    @property
    def start(self):
        return min(e.start for e in self.exons)

    @property
    def stop(self):
        return max(e.stop for e in self.exons)

    @property
    def string_strand(self):
        if self.strand == 1:
            return '+'
        if self.strand == -1:
            return '-'
        raise UnknownStrand()

    @property
    def name(self):
        exon_names = []
        for exon in self.exons:
            exon_names.append('{start}-{stop}'.format(
                start=exon.start,
                stop=exon.stop,
            ))
        return '@{chromosome}/{exons}:{strand}'.format(
            chromosome=self.chromosome,
            exons=','.join(exon_names),
            strand=self.string_strand,
        )

    def writeable(self, accession):
        for exon in self.exons:
            yield [
                accession,
                self.name,
                self.chromosome,
                self.strand,
                self.assembly_id,
                len(self.exons),
                exon.start,
                exon.stop,
            ]

    def writeable_exons(self, accession):
        for exon in self.exons:
            yield [
                accession,
                self.chromosome,
                exon.start,
                exon.stop,
                self.assembly_id,
                self.strand,
            ]
