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

import re
import operator as op
import collections as coll

from rnacentral_pipeline.databases.helpers.publications import reference
from ..data import Exon

seq_version = op.itemgetter('seq_version')
species = op.itemgetter('species')
product = op.itemgetter('description')
primary_id = op.itemgetter('primary_id')
common_name = op.itemgetter('common_name')
optional_id = op.itemgetter('optional_id')
parent_accession = op.itemgetter('parent_accession')


def lineage(data):
    base = re.sub(r'\.$', '', data['lineage'])
    return '; '.join(base.split(' ') + [data['species']])


def taxid(data):
    return int(data['ncbi_tax_id'])


def feature_location_endpoints(data):
    start = int(data['feature_location_start'])
    stop = int(data['feature_location_end'])
    return (start, stop)



def mol_type(data):
    if data['is_seed'] == '0':
        return 'full'
    return 'seed'


def references(data):
    refs = [reference(29112718)]
    refs.extend(reference(pmid) for pmid in data.get('references', []))
    return refs


def note(data):
    result = coll.defaultdict(list)
    result['Alignment'] = mol_type(data)
    for entry in data['ontology']:
        db = entry.split(':', 1)[0]
        result[db].append(entry)
    return dict(result)


def exons(data):
    location_range = feature_location_endpoints(data)
    complement = int(data['feature_location_start']) != location_range[0]
    return [Exon(
        chromosome_name='',
        primary_start=location_range[0],
        primary_end=location_range[1],
        assembly_id='',
        complement=complement,
    )]


def rna_type(data, mapping):
    return mapping.get(data['primary_id'], None) or data['ncrna_class']


def sequence(data):
    return data['sequence'].upper()


def experiment(data):
    return ' '.join(data['references'])


def accession(data):
    location_range = feature_location_endpoints(data)
    return ('{parent}.{version}:{start}..{stop}:rfam').format(
        parent=parent_accession(data),
        version=seq_version(data),
        start=location_range[0],
        stop=location_range[1],
    )


def description(data):
    return '{species} {product}'.format(
        species=species(data),
        product=product(data)
    )


def url(data):
    return 'http://rfam.org/family/%s' % data['primary_id']
