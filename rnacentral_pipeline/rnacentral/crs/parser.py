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
import json
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a


@attr.s(frozen=True)
class GenomicLocation(object):
    chromosome = attr.ib(validator=is_a(basestring))
    start = attr.ib(validator=is_a(int), convert=int)
    stop = attr.ib(validator=is_a(int), convert=int)

    @classmethod
    def build(cls, raw):
        return cls(
            chromosome=raw['chromosome'],
            start=raw['start'],
            stop=raw['stop'],
        )


@attr.s(frozen=True)
class CrsFeature(object):
    upi = attr.ib(validator=is_a(basestring))
    taxid = attr.ib(validator=is_a(int), convert=int)
    crs_name = attr.ib(validator=is_a(basestring))
    start = attr.ib(validator=is_a(int), convert=int)
    stop = attr.ib(validator=is_a(int), convert=int)
    genomic_locations = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, raw_features):
        upi, taxid = raw_features[0]['URS_taxid'].split('_')
        return cls(
            upi=upi,
            taxid=taxid,
            crs_name=raw_features[0]['CRS_id'],
            start=raw_features[0]['start'],
            stop=raw_features[0]['stop'],
            genomic_locations=[GenomicLocation.build(f) for f in raw_features],
        )

    def writeable(self):
        metadata = attr.asdict(self)
        metadata = {'genomic_locations': metadata['genomic_locations']}
        return [
            self.upi,
            self.taxid,
            None,
            self.start,
            self.stop,
            'conserved_rna_structure',
            json.dumps(metadata),
        ]


def parse(handle):
    key = op.itemgetter('URS_taxid', 'CRS_id', 'start', 'stop')
    raw = sorted(csv.DictReader(handle, delimiter='\t'), key=key)
    for _, features in it.groupby(raw, key):
        yield CrsFeature.build(list(features))


def from_file(handle, output):
    data = parse(handle)
    data = it.imap(op.methodcaller('writeable'), data)
    writer = csv.writer(output, delimiter=',')
    writer.writerows(data)
