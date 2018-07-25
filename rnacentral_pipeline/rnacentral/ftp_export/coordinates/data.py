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

import json
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a
from attr.validators import in_ as one_of


@attr.s()
class ExonicSequence(object):
    rna_id = attr.ib(validator=is_a(basestring))
    rna_type = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    exons = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, raw):
        exons = []
        if raw['known_coordinates']['region_id'] and \
                raw['known_coordinates']['start'] is not None:
            exons.append(Exon.build(raw['known_coordinates']))
        if raw['mapped_coordinates']['region_id'] and \
                raw['mapped_coordinates']['start'] is not None:
            exons.append(Exon.build(raw['mapped_coordinates']))

        assert exons, "Could build any exons"
        return cls(
            rna_id=raw['rna_id'],
            rna_type=raw['rna_type'],
            databases=raw['databases'],
            exons=exons,
        )


@attr.s(slots=True, frozen=True)
class LocatedSequence(object):
    rna_id = attr.ib(validator=is_a(basestring))
    rna_type = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(basestring))
    regions = attr.ib(validator=is_a(list))

    @classmethod
    def from_exonics(cls, exonics):
        exons = []
        exonics = list(exonics)
        for exonic in exonics:
            exons.extend(exonic.exons)

        regions = []
        key = op.attrgetter('region_id')
        exon_key = op.attrgetter('chromosome', 'start', 'stop')
        exons.sort(key=key)
        seen = set()
        for index, (_, exons) in enumerate(it.groupby(exons, key)):
            ordered = sorted(exons, key=exon_key)
            region_id = '%s.%i' % (exonics[0].rna_id, index)
            region = Region.from_exons(region_id, ordered)
            if region not in seen:
                regions.append(region)
                seen.add(region)

        return cls(
            rna_id=exonics[0].rna_id,
            rna_type=exonics[0].rna_type,
            databases=exonics[0].databases,
            regions=regions,
        )


@attr.s(hash=True, slots=True, frozen=True)
class Region(object):
    region_id = attr.ib(validator=is_a(basestring), hash=False, cmp=False)
    chromosome = attr.ib(validator=is_a(basestring))
    strand = attr.ib(validator=is_a(int))
    endpoints = attr.ib(validator=is_a(tuple))
    coordinate_system = attr.ib(validator=one_of(['chromosome', 'scaffold']))

    @classmethod
    def from_exons(cls, region_id, exons):
        region_ids = set()
        chromosomes = set()
        strands = set()
        endpoints = set()

        for exon in exons:
            region_ids.add(exon.region_id)
            chromosomes.add(exon.chromosome)
            strands.add(exon.strand)
            endpoints.add(Endpoint.from_exon(exon))

        assert len(region_ids) == 1
        assert len(chromosomes) == 1
        assert len(strands) == 1

        return cls(
            region_id=region_id,
            chromosome=chromosomes.pop(),
            strand=strands.pop(),
            endpoints=tuple(sorted(endpoints, key=op.attrgetter('start', 'stop')))
        )

    @property
    def start(self):
        return self.endpoints[0].start

    @property
    def stop(self):
        return self.endpoints[-1].stop

    def string_strand(self):
        if self.strand == 1:
            return '+'
        if self.strand == -1:
            return '-'
        raise ValueError("Unknown type of strand")

    def is_chromosomal(self):
        return self.coordinate_system == 'chromosome'


@attr.s()
class Exon(object):
    region_id = attr.ib(validator=is_a(basestring))
    chromosome = attr.ib(validator=is_a(basestring))
    strand = attr.ib(validator=is_a(int), convert=int)
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, raw):
        return cls(**raw)


@attr.s(hash=True, slots=True, frozen=True)
class Endpoint(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def from_exon(cls, exon):
        return cls(start=exon.start, stop=exon.stop)


def is_buildable(raw):
    return (
        raw['known_coordinates']['region_id'] and
        raw['known_coordinates']['start'] is not None
    ) or (
        raw['mapped_coordinates']['region_id'] and
        raw['mapped_coordinates']['start'] is not None
    )

def parse(iterable):
    buildable = it.ifilter(is_buildable, iterable)
    exonics = it.imap(ExonicSequence.build, buildable)
    grouped = it.groupby(exonics, op.attrgetter('rna_id'))
    grouped = it.imap(op.itemgetter(1), grouped)
    located = it.imap(LocatedSequence.from_exonics, grouped)
    return located


def from_file(handle):
    data = it.imap(json.loads, handle)
    return parse(data)
