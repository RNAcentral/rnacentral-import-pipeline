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

import collections as coll
import re

import requests
from functools32 import lru_cache

TAX_URL = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{taxon_id}'


CODING_RNA_TYPES = set([
    'tec',
    'processed_transcript',
    'retained_intron',
])


class MissingTaxId(Exception):
    pass


def qualifier_value(feature, name, pattern, max_allowed=1):
    values = set()
    for note in feature.qualifiers.get(name, []):
        match = re.match(pattern, note)
        if match:
            values.add(match.group(1))
    if max_allowed is not None and len(values) > max_allowed:
        raise ValueError("Multiple values (%s) for %s",
                         ', '.join(sorted(values)), name)

    if len(values) == 0:
        return None
    if max_allowed == 1:
        return values.pop()
    return values


def taxid(record):
    source = record.features[0]
    if source.type == 'source':
        return int(qualifier_value(source, 'db_xref', r'^taxon:(\d+)$'))
    raise MissingTaxId("No taxon id found for record %s" % record)


def organism_naming(record):
    pattern = re.compile(r'\((.+)\)$')
    common_name = None
    species = record.annotations['organism']
    match = re.search(pattern, species)
    if match:
        common_name = match.group(1)
        species = re.sub(pattern, '', species).strip()
    return (species, common_name)


def gene(feature):
    return qualifier_value(feature, 'gene', '^(.+)$')


def is_gene(feature):
    return feature.type == 'gene'


def transcript(feature):
    return qualifier_value(feature, 'note', '^transcript_id=(.+)$')


def is_transcript(feature):
    return transcript(feature) is not None


def notes(feature):
    return feature.qualifiers.get('note', [])


def rna_type(feature):
    return notes(feature)[0]


def locus_tag(feature):
    return feature.qualifiers.get('locus_tag', [''])[0]


def note_data(feature):
    data = notes(feature)
    if len(data) > 1:
        data = data[1:]
    return grouped_annotations(data, '=')


def xref_data(feature):
    raw = feature.qualifiers.get('db_xref', [])
    return grouped_annotations(raw, ':')


def grouped_annotations(raw, split):
    parsed = coll.defaultdict(set)
    for entry in raw:
        if split not in entry:
            continue
        key, value = entry.split(split, 1)
        if key == 'RNACentral':
            key = 'RNAcentral'
        parsed[key].add(value)
    return {k: sorted(v) for k, v in parsed.items()}


def is_ncrna(feature):
    """
    Checks if the feature is that of ncRNA. This requires that the type be
    one of misc_RNA or ncRNA and that the ncRNA type is not 'TEC'.
    """

    # The checking the first entry in 'note' is a quick and dirty way of
    # getting the ncRNA type.
    return feature.type in {'misc_RNA', 'ncRNA'} and \
        feature.qualifiers['note'][0].lower() not in CODING_RNA_TYPES


def standard_name(feature):
    return qualifier_value(feature, 'standard_name', '^(.+)$')


def chromosome(record):
    basic = record.id
    if basic.startswith('chromosome:') or basic.startswith('scaffold:'):
        parts = basic.split(':')
        return parts[2]
    return basic.split('.')[0]


@lru_cache()
def lineage(record):
    taxon_id = taxid(record)
    if taxon_id <= 0:
        return None

    response = requests.get(TAX_URL.format(taxon_id=taxon_id))
    response.raise_for_status()
    response_data = response.json()
    return '{lineage}{name}'.format(
        lineage=response_data['lineage'],
        name=response_data['scientificName']
    )
