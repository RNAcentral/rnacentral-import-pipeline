# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

from jsonschema import validate


MOD_URL = 'http://modomics.genesilico.pl/sequences/list/{id}'


TRUSTED_DB = set([
    'gtrnadb',
    'lncrnadb',
    'mirbase',
    'modomics',
    'pdbe',
    'snopy'
    'srpdb',
    'tmrna website',
])


def external_id(data):
    if data['database'] == 'PDBe':
        return '%s_%s' % (data['external_id'], data['optional_id'])
    if data['database'] == 'Modomics':
        return MOD_URL.format(id=data['external_id'])
    return data['external_id']


def is_high_quality(data):
    name = data['database'].lower()
    if name in TRUSTED_DB:
        return True
    if name == 'rfam':
        return data['molecule_type'] == 'seed'
    return False


def as_xref(xref):
    return {'database': xref['database'], 'id': external_id(xref)}


def builder(data):
    result = dict(data)
    xrefs = []
    seen = set()
    key = op.itemgetter('database', 'id')
    for xref in data['xrefs']:
        if not is_high_quality(xref):
            continue

        updated = as_xref(xref)
        value = key(updated)
        if value not in seen:
            xrefs.append(updated)
            seen.add(value)

    result['sequence'] = result['sequence'].upper().replace('U', 'T')
    result['xrefs'] = xrefs
    return result


def generate_file(raw, output, schema_file=None):
    results = [builder(json.loads(l)) for l in raw]

    if schema_file:
        with open(schema_file, 'r') as raw:
            validate(results, json.load(raw))

    json.dump(results, output)
