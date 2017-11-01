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
import collections as coll

import requests
from functools32 import lru_cache

from databases.data import Exon

TAX_URL = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/tax-id/{taxon_id}'


CODING_RNA_TYPES = set([
    'tec',
    'processed_transcript',
    'retained_intron',
])


class MissingTaxId(Exception):
    """
    This is  raised when an operation which should have a NCBI taxon id should,
    but does not.
    """
    pass


def qualifier_value(feature, name, pattern, max_allowed=1):
    """
    This will parse the qualifer feild defined by the given name for the given
    feature. This will extract all values matching the given regex pattern. If
    max allowed is 1 then only one distinct value is allowed and a single value
    will be returned. Otherwise all values in a set will be returned.
    """

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
    """
    Get the taxon id of the given record. This will pull the first feature,
    which must be of the 'source' type to do so.
    """

    source = record.features[0]
    if source.type == 'source':
        return int(qualifier_value(source, 'db_xref', r'^taxon:(\d+)$'))
    raise MissingTaxId("No taxon id found for record %s" % record)


def organism_naming(record):
    """
    This will parse out the species and common name from the assigned species
    of this record. The format is 'Homo sapiens (human)' where the string in
    the parenthesis is the common name and everything else is the species.
    """

    pattern = re.compile(r'\((.+)\)$')
    common_name = None
    species = record.annotations['organism']
    match = re.search(pattern, species)
    if match:
        common_name = match.group(1)
        species = re.sub(pattern, '', species).strip()
    return (species, common_name)


def gene(feature):
    """
    Get the gene this feature is a part of.
    """
    return qualifier_value(feature, 'gene', '^(.+)$')


def is_gene(feature):
    """
    Check if this feature is a gene
    """
    return feature.type == 'gene'


def transcript(feature):
    """
    Get the transcript id of this feature.
    """
    return qualifier_value(feature, 'note', '^transcript_id=(.+)$')


def is_transcript(feature):
    """
    Check if this feature is a transcript.
    """
    return transcript(feature) is not None


def notes(feature):
    """
    Get all notes from this feature. If none a present then an empty list is
    returned.
    """
    return feature.qualifiers.get('note', [])


def rna_type(feature):
    """
    Get the RNA type from the feature.
    """
    return notes(feature)[0]


def locus_tag(feature):
    """
    Get the locus tag of this feature. If none is present then the empty string
    is returned.
    """
    return feature.qualifiers.get('locus_tag', [''])[0]


def note_data(feature):
    """
    This will parse the notes data of the feature to produce a dict of key
    value mappings. This will strip out the RNA type in the final dictonary.
    """
    notes_data = notes(feature)
    if len(notes_data) > 1:
        notes_data = notes_data[1:]
    return grouped_annotations(notes_data, '=')


def xref_data(feature):
    """
    Get a dict of the xref data for this feature. This will parse the db_xref
    qualifier to produce a key value mapping.
    """
    raw = feature.qualifiers.get('db_xref', [])
    return grouped_annotations(raw, ':')


def grouped_annotations(raw, split):
    """
    Parse a raw string into a dict. This will produce a key value mappign where
    the key is everything before the first split and values is everything
    after. The mapping values will be lists. This will correct RNACentral to
    RNAcentral.
    """
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
    """
    Get the standard name of feature.
    """
    return qualifier_value(feature, 'standard_name', '^(.+)$')


def chromosome(record):
    """
    Get the chromosome this record is from. If it is part of a scaffold then it
    will be returned.
    """
    basic = record.id
    if basic.startswith('chromosome:') or basic.startswith('scaffold:'):
        parts = basic.split(':')
        return parts[2]
    return basic.split('.')[0]


@lru_cache()
def lineage(record):
    """
    Extract the taxon id and then query a remote server for the lineage for the
    given taxon id. Results are cached because it is common to constantly query
    with only the same few taxon ids.
    """
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


def product(feature):
    """
    Determine the product of this feature. If this has RNA type of scaRNA
    then we use that as the product, otherwise nothing.
    """

    if rna_type(feature) == 'scaRNA':
        return 'scaRNA'
    return None


def exon(location):
    """
    Build an Exon from a biopython location object.
    """
    return Exon(
        chromosome='',
        primary_start=location.start + 1,
        primary_end=int(location.end),
        complement=location.strand == -1,
    )
