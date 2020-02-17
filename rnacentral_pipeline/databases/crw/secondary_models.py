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
import operator as op

import six
import attr
from attr.validators import instance_of as is_a

import typing

SO_TERM_MAPPING = {
    '16S': 'SO:0000650',
    '23S': 'SO:0000651',
    '5S': 'SO:0000652',
    'I': 'SO:0000587',
    'IA1': 'SO:0000587',
    'IA2': 'SO:0000587',
    'IB': 'SO:0000587',
    'IB1': 'SO:0000587',
    'IB2': 'SO:0000587',
    'IB4': 'SO:0000587',
    'IC1': 'SO:0000587',
    'IC2': 'SO:0000587',
    'IC3': 'SO:0000587',
    'ID': 'SO:0000587',
    'IE': 'SO:0000587',
    'IIA': 'SO:0000603',
    'IIB': 'SO:0000603',
}


def as_so_term(raw):
    if raw in SO_TERM_MAPPING:
        return SO_TERM_MAPPING[raw]
    raise ValueError("Unknown RNA type: " + raw)


def as_taxid(raw):
    if raw == '501083':
        return 126
    if raw in {'600001', '600002', '600003'}:
        return 562
    if raw in {'600101', '600102'}:
        return 2238
    if raw in {'600301', '600302'}:
        return 4932
    if raw in {'600201', '600202'}:
        return 274
    return int(raw)


@attr.s()
class Info(object):
    model_id = attr.ib(validator=is_a(six.text_type))
    is_intronic = attr.ib(validator=is_a(bool))
    so_term = attr.ib(validator=is_a(six.text_type))
    taxid = attr.ib(validator=is_a(six.integer_types))
    accessions = attr.ib(type=typing.List[six.text_type])
    cell_location = attr.ib(validator=is_a(six.text_type))

    @classmethod
    def build(cls, model_id, raw):
        intronic = raw['rna_type'] == 'I'
        return cls(
            model_id=model_id,
            is_intronic=intronic,
            so_term=as_so_term(raw['rna_class']),
            taxid=as_taxid(raw['tax_id']),
            accessions=raw['accession(s)'].split(','),
            cell_location=raw['cell_location'],
        )

    @classmethod
    def build_all(cls, raw):
        for model_id in raw['structure'].split(' '):
            model_id = re.sub(r'\.ps$', '', model_id)
            yield cls.build(model_id, raw)

    @property
    def rna_type(self):
        if self.so_term in {'SO:0000650', 'SO:0000651', 'SO:0000652'}:
            return 'rRNA'
        if self.so_term in {'SO:0000587', 'SO:0000603'}:
            return 'autocatalytically_spliced_intron'
        raise ValueError("No RNA type for: " + self.so_term)

    def writeable(self):
        return [
            self.model_id,
            self.taxid,
            self.rna_type,
            self.so_term,
            self.cell_location,
        ]


def parse(handle):
    for row in csv.DictReader(handle, delimiter='\t'):
        for info in Info.build_all(row):
            yield info


def write(handle, output):
    data = parse(handle)
    data = six.moves.map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
