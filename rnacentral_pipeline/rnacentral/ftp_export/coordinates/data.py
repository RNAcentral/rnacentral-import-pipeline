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
from attr.validators import optional
from attr.validators import in_ as one_of
from attr.validators import instance_of as is_a


@attr.s()
class ExonicSequence(object):
    rna_id = attr.ib(validator=is_a(basestring))
    rna_type = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(list))
    exons = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, raw):
        exons = [Exon.build(e) for e in raw['exons']]
        assert exons, "Could build any exons"

        if raw['rna_type'] in {'miRNA', 'piRNA'}:
            updates = []
            for index, exon in enumerate(exons):
                region_id = exon.region_id + ':%s' % (index + 1)
                updates.append(attr.evolve(exon, region_id=region_id))
            exons = updates
        exons = sorted(set(exons))

        return cls(
            rna_id=raw['rna_id'],
            rna_type=raw['rna_type'],
            databases=raw['databases'].split(','),
            exons=exons,
        )


@attr.s(slots=True, frozen=True)
class LocatedSequence(object):
    rna_id = attr.ib(validator=is_a(basestring))
    rna_type = attr.ib(validator=is_a(basestring))
    databases = attr.ib(validator=is_a(list))
    regions = attr.ib(validator=is_a(list))

    @classmethod
    def from_exonic(cls, exonic):
        seen = set()
        regions = []
        grouping_key = op.attrgetter('region_id')
        exons = list(exonic.exons)
        exons.sort(key=grouping_key)
        for _, exons in it.groupby(exons, grouping_key):
            ordered = sorted(exons, key=op.attrgetter('start', 'stop'))

            sources = {e.source for e in ordered}
            if not sources:
                raise ValueError("No sources defined")
            elif len(sources) == 1:
                source = sources.pop()
            elif sources == {'expert-database', 'alignment'}:
                source = 'expert-database'
            else:
                raise ValueError("Unknown sources in: %s" % str(sources))

            region = Region.from_exons(
                exonic.rna_id,
                source,
                ordered
            )
            if region not in seen:
                regions.append(region)
                seen.add(region)

        def region_key(region):
            try:
                chromosome = int(region.chromosome)
            except:
                chromosome = region.chromosome
            return (chromosome, region.strand, region.start, region.stop)

        regions.sort(key=region_key)
        return cls(
            rna_id=exonic.rna_id,
            rna_type=exonic.rna_type,
            databases=exonic.databases,
            regions=regions,
        )


@attr.s(hash=True, slots=True, frozen=True)
class Region(object):
    rna_id = attr.ib(validator=is_a(basestring))
    chromosome = attr.ib(validator=is_a(basestring))
    strand = attr.ib(validator=is_a(int))
    endpoints = attr.ib(validator=is_a(tuple))
    source = attr.ib(
        validator=one_of(['expert-database', 'alignment']),
        cmp=False,
    )
    # coordinate_system = attr.ib(validator=one_of(['chromosome', 'scaffold']))
    identity = attr.ib(
        validator=optional(is_a(float)),
        default=None,
        cmp=False,
    )

    @classmethod
    def from_exons(cls, rna_id, source, exons):
        region_ids = set()
        chromosomes = set()
        strands = set()
        endpoints = set()

        for exon in exons:
            region_ids.add(exon.region_id)
            chromosomes.add(exon.chromosome)
            strands.add(exon.strand)
            endpoints.add(Endpoint.from_exon(exon))

        if len(chromosomes) != 1:
            from pprint import pprint
            pprint(rna_id)
            pprint(source)
            pprint(exons)

        assert len(region_ids) == 1
        assert len(chromosomes) == 1
        assert len(strands) == 1

        identity = exons[0].identity
        if source == 'expert-database':
            identity = None

        return cls(
            rna_id=rna_id,
            chromosome=chromosomes.pop(),
            strand=strands.pop(),
            source=source,
            endpoints=tuple(sorted(endpoints, key=op.attrgetter('start', 'stop'))),
            identity=identity
        )

    @property
    def region_id(self):
        return '{rna_id}@{chromosome}/{start}-{stop}:{strand}'.format(
            rna_id=self.rna_id,
            chromosome=self.chromosome,
            start=self.start,
            stop=self.stop,
            strand=self.string_strand(),
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

    # def is_chromosomal(self):
    #     return self.coordinate_system == 'chromosome'


@attr.s(hash=True)
class Exon(object):
    region_id = attr.ib(validator=is_a(basestring))
    chromosome = attr.ib(validator=is_a(basestring))
    strand = attr.ib(validator=is_a(int), convert=int)
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))
    source = attr.ib(validator=is_a(basestring))
    identity = attr.ib(validator=optional(is_a(float)))

    @classmethod
    def build(cls, raw):
        identity = None
        if 'identity' in raw:
            identity = float(raw['identity'])

        region_id = raw['region_id']
        if raw['source'] == 'expert-database':
            region_id = '{region_id}:{chr}:{strand}'.format(
                region_id=region_id,
                chr=raw['chromosome'],
                strand=raw['strand'],
            )

        return cls(
            region_id=region_id,
            chromosome=raw['chromosome'],
            strand=raw['strand'],
            start=raw['start'],
            stop=raw['stop'],
            source=raw['source'],
            identity=identity
        )


@attr.s(hash=True, slots=True, frozen=True)
class Endpoint(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def from_exon(cls, exon):
        return cls(start=exon.start, stop=exon.stop)


def is_buildable(raw):
    return raw['region_id'] and raw['start'] is not None

def parse(iterable):
    sources = {
        'known_coordinates': 'expert-database',
        'mapped_coordinates': 'alignment',
    }

    for _, entries in it.groupby(iterable, op.itemgetter('rna_id')):
        exons = []
        entries = list(entries)
        for entry in entries:
            common = {k: v for k, v in entry.items() if k not in sources}
            for name, source in sources.items():
                data = entry[name]
                if not is_buildable(data):
                    continue
                current = {'source': source}
                current.update(common)
                current.update(data)
                exons.append(current)

        exonic = dict(common)
        exonic['exons'] = exons
        exonic = ExonicSequence.build(exonic)
        yield LocatedSequence.from_exonic(exonic)


def from_file(handle):
    data = it.imap(json.loads, handle)
    return parse(data)
