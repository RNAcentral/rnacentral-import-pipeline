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

from rnacentral_pipeline.rnacentral.traveler.data import ModelInfo

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


def models(raw):
    for model_id in raw['structure'].split(' '):
        data = dict(raw)
        model_id = re.sub(r'\.ps$', '', model_id)
        data['model_id'] = model_id
        yield data


def parse(handle):
    for row in csv.DictReader(handle, delimiter='\t'):
        for info in models(row):
            intronic = info['rna_type'] == 'I'
            yield ModelInfo(
                model_id=info['model_id'],
                is_intronic=intronic,
                so_term=as_so_term(info['rna_class']),
                taxid=as_taxid(info['tax_id']),
                accessions=row['accession(s)'].split(','),
                cell_location=info['cell_location'],
            )


def write(handle, output):
    data = parse(handle)
    data = map(op.methodcaller('writeable'), data)
    csv.writer(output).writerows(data)
