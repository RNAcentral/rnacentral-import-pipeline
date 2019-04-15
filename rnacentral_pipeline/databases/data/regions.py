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

import operator as op

import six

import enum

import attr
from attr.validators import instance_of as is_a


class UnknownStrand(Exception):
    """
    Raised when a strand integer has an invalid value.
    """
    pass


def as_strand(value):
    if isinstance(value, six.integer_types):
        return value
    elif isinstance(value, float) and int(value) == value:
        return int(value)
    else:
        if value == '+' or value == '1':
            return 1
        elif value == '-' or value == '-1':
            return -1
        elif value == '.':
            return 0
    raise UnknownStrand("No way to handle strand: " + str(value))


def sort_exons(exons):
    return sorted(exons, key=op.attrgetter('start'))


def string_strand(located):
    if located.strand == 1 or located.strand == '+':
        return '+'
    if located.strand == -1 or located.strand == '-':
        return '-'
    if located.strand == 0 or located.strand == '.':
        return '.'
    raise UnknownStrand()


def region_id(located):
    exon_names = []
    for exon in located.exons:
        exon_names.append('{start}-{stop}'.format(
            start=exon.start,
            stop=exon.stop,
        ))
    return '{upi}@{chromosome}/{exons}:{strand}'.format(
        upi=getattr(located, 'upi', ''),
        chromosome=located.chromosome,
        exons=','.join(exon_names),
        strand=string_strand(located),
    )


def write_locations(located, accession, require_strand=True):
    if require_strand and located.strand == 0:
        return

    for exon in located.exons:
        normalized = as_1_based(exon, located.coordinate_system)
        yield [
            accession,
            located.name,
            located.chromosome,
            located.strand,
            located.assembly_id,
            len(located.exons),
            normalized.start,
            normalized.stop,
        ]


def as_1_based(coordinate_system, exon):
    if coordinate_system == CoordinateSystem.one_based:
        return exon
    elif coordinate_system == CoordinateSystem.zero_based:
        return None
    raise ValueError("Unknown coordinate_system: %s" % coordinate_system)


@enum.unique
class CoordinateSystem(enum.Enum):
    zero_based = 0
    one_based = 1


@attr.s(frozen=True)
class Exon(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))


@attr.s()
class SequenceRegion(object):
    assembly_id = attr.ib(validator=is_a(six.text_type))
    chromosome = attr.ib(validator=is_a(six.text_type))
    strand = attr.ib(validator=is_a(int), converter=as_strand)
    exons = attr.ib(validator=is_a(list), converter=sort_exons)
    coordinate_system = attr.ib(validator=is_a(CoordinateSystem))

    @property
    def start(self):
        return min(e.start for e in self.exons)

    @property
    def stop(self):
        return max(e.stop for e in self.exons)

    @property
    def name(self):
        return region_id(self)

    def writeable(self, accession):
        return write_locations(self, accession)
