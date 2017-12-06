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

from databases.data import Exon
from databases.helpers.embl import taxid
from databases.helpers.embl import lineage
from databases.helpers.embl import qualifier_value
from databases.helpers.embl import standard_name
from databases.helpers.embl import gene
from databases.helpers.embl import xref_data
from databases.helpers.embl import grouped_annotations


CODING_RNA_TYPES = set([
    'tec',
    'processed_transcript',
    'retained_intron',
])


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


def note_data(feature):
    """
    This will parse the notes data of the feature to produce a dict of key
    value mappings. This will strip out the RNA type in the final dictonary.
    """
    notes_data = notes(feature)
    if len(notes_data) > 1:
        notes_data = notes_data[1:]
    return grouped_annotations(notes_data, '=')


def is_ncrna(feature):
    """
    Checks if the feature is that of ncRNA. This requires that the type be
    one of misc_RNA or ncRNA and that the ncRNA type is not 'TEC'.
    """

    # The checking the first entry in 'note' is a quick and dirty way of
    # getting the ncRNA type.
    return feature.type in {'misc_RNA', 'ncRNA'} and \
        feature.qualifiers['note'][0].lower() not in CODING_RNA_TYPES


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
