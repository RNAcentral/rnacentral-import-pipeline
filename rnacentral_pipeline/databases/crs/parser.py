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

import re
import csv
import json
import operator as op
import itertools as it

import attr
from attr.validators import instance_of as is_a

import six


@attr.s(frozen=True)
class GenomicLocation(object):
    chromosome = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int), converter=int)
    stop = attr.ib(validator=is_a(int), converter=int)

    @classmethod
    def build(cls, raw):
        chromosome = re.sub('^chr', '', raw['chromosome'])
        return cls(
            chromosome=chromosome,
            start=raw['CRS_start_relative_to_genome'],
            stop=raw['CRS_end_relative_to_genome'],
        )


@attr.s(frozen=True)
class CrsFeature(object):
    upi = attr.ib(validator=is_a(str))
    taxid = attr.ib(validator=is_a(int), converter=int)
    crs_name = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int), converter=int)
    stop = attr.ib(validator=is_a(int), converter=int)
    genomic_location = attr.ib(validator=is_a(GenomicLocation))
    fdr = attr.ib(validator=is_a(float), converter=float)

    @classmethod
    def build(cls, raw_feature):
        upi, taxid = raw_feature['URS_taxid'].split('_')
        return cls(
            upi=upi,
            taxid=taxid,
            crs_name=raw_feature['CRS_id'],
            start=raw_feature['CRS_start_relative_to_URS'],
            stop=raw_feature['CRS_end_relative_to_URS'],
            genomic_location=GenomicLocation.build(raw_feature),
            fdr=raw_feature['CRS_fdr'],
        )

    def writeable(self):
        metadata = attr.asdict(self)
        metadata = {
            'crs_id': self.crs_name,
            'genomic_location': metadata['genomic_location'],
            'fdr': self.fdr,
            'should_highlight': self.fdr <= 10.0,
        }

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
    field_names = [
        'URS_taxid',
        'CRS_start_relative_to_URS',
        'CRS_end_relative_to_URS',
        'chromosome',
        'CRS_start_relative_to_genome',
        'CRS_end_relative_to_genome',
        'CRS_id',
        'CRS_fdr',
    ]
    for row in csv.reader(handle, delimiter='\t'):
        if row[0] == 'URS_taxid':
            continue
        assert len(row) == len(field_names)
        feature = dict(zip(field_names, row))
        yield CrsFeature.build(feature)


def from_file(handle, output):
    data = parse(handle)
    data = six.moves.map(op.methodcaller('writeable'), data)
    writer = csv.writer(output, delimiter=',')
    writer.writerows(data)
