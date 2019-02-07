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
from attr.validators import instance_of as is_a


def clean_databases(raw):
    return [d.replace(' ', '_') for d in raw]


@attr.s(hash=True, slots=True, frozen=True)
class Region(object):
    region_id = attr.ib(validator=is_a(str))
    rna_id = attr.ib(validator=is_a(str))
    chromosome = attr.ib(validator=is_a(str))
    strand = attr.ib(validator=is_a(int))
    endpoints = attr.ib(validator=is_a(tuple))
    was_mapped = attr.ib(validator=is_a(bool))
    identity = attr.ib(
        validator=optional(is_a(float)),
        default=None,
        cmp=False,
    )
    metadata = attr.ib(validator=is_a(dict), default=dict, cmp=False)

    @classmethod
    def build(cls, index, raw):
        identity = None
        if raw['identity'] is not None:
            identity = float(raw['identity'])

        endpoints = [Endpoint.build(e) for e in raw['exons']]
        endpoints.sort(key=op.attrgetter('start'))
        endpoints = tuple(endpoints)

        region_id = '{rna_id}.{index}'.format(
            rna_id=raw['rna_id'],
            index=index
        )

        metadata = {
            'rna_type': raw['rna_type'],
            'providing_databases': clean_databases(raw['providing_databases']),
            'databases': clean_databases(raw['databases']),
        }
        if not metadata['providing_databases']:
            del metadata['providing_databases']

        return cls(
            region_id=region_id,
            rna_id=raw['rna_id'],
            chromosome=raw['chromosome'],
            strand=raw['strand'],
            endpoints=endpoints,
            identity=identity,
            was_mapped=raw['was_mapped'],
            metadata=metadata,
        )

    @property
    def start(self):
        return self.endpoints[0].start

    @property
    def stop(self):
        return self.endpoints[-1].stop

    @property
    def source(self):
        if self.was_mapped:
            return 'alignment'
        return 'expert-database'

    def string_strand(self):
        if self.strand == 1:
            return '+'
        if self.strand == -1:
            return '-'
        raise ValueError("Unknown type of strand")


@attr.s(hash=True, slots=True, frozen=True)
class Endpoint(object):
    start = attr.ib(validator=is_a(int))
    stop = attr.ib(validator=is_a(int))

    @classmethod
    def build(cls, raw):
        return cls(start=raw['exon_start'], stop=raw['exon_stop'])


def parse(regions):
    for index, region in enumerate(regions):
        yield Region.build(index, region)


def from_file(handle):
    data = it.imap(json.loads, handle)
    return parse(data)
