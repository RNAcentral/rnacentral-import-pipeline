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
import typing as ty

import attr

from rnacentral_pipeline.databases import data

from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.helpers import publications as pubs

from rnacentral_pipeline.databases.ensembl.vertebrates.context import Context

URL = 'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript}'

CODING_RNA_TYPES = set([
    'tec',
    'processed_transcript',
    'retained_intron',
    'crispr',
])


class CouldNotTrimDescription(Exception):
    """
    Raised when a description could not be trimmed to only the organisms
    suffix.
    """
    pass


def references():
    """
    Get the standard reference for all Ensembl entries.
    """
    return [pubs.reference(27337980)]


def is_pseudogene(gene, feature):
    """
    Check if the given feature is a pseudogene. This is determined by check
    if the feature is labeled as a psuedogene or the gene it is a transcript
    from a pseudogene.
    """

    raw_notes = feature.qualifiers.get('note', [])
    raw_notes.extend(gene.qualifiers.get('note', []))
    return any('pseudogene' in n for n in raw_notes)


def raw_rna_type(feature):
    """
    Get the RNA type from the feature.
    """
    return notes(feature)[0]


def rna_type(inference, feature, xref_data):
    """
    Compute the RNA type of the given feature.
    """

    base_type = raw_rna_type(feature)
    return inference.infer_rna_type(xref_data, base_type)


def description(context, gene, entry):
    """
    Generate a description for the entry based upon the locus this is a part
    of. This will be of the form 'Homo sapiens (human) lncRNA xist
    """

    species = entry.species
    if entry.common_name:
        species += ' (%s)' % entry.common_name

    gene_name = notes(gene)
    if gene_name and len(gene_name) == 1:
        gene_name = gene_name[0]
        if gene_name.endswith(']'):
            gene_name = re.sub(r'\s*\[.+\]$', '', gene_name)
        gene_name.strip()
        return '{species} {gene_name}'.format(
            species=species,
            gene_name=gene_name
        )

    locus_tag = context.rfam_name(entry.locus_tag, entry.locus_tag or '')
    if locus_tag.startswith('RF') and entry.optional_id:
        locus_tag = entry.optional_id.split('.')[0]

    assert entry.rna_type, "Cannot build description without rna_type"
    rna_type = entry.human_rna_type()
    if rna_type == locus_tag:
        locus_tag = ''
    return '{species} {rna_type} {locus_tag}'.format(
        species=species,
        rna_type=rna_type,
        locus_tag=locus_tag,
    ).strip()


def seq_version(feature):
    """
    Compute the sequence version, if any, of the given id.
    """

    acc = accession(feature)
    if '.' in acc:
        return acc.split('.', 1)[1]
    return '1'


def transcript(feature):
    """
    Get the transcript id of this feature.
    """
    name = embl.qualifier_value(feature, 'note', '^transcript_id=(.+)$')
    if name:
        return name
    return embl.qualifier_string(feature, 'standard_name')


def primary_id(feature):
    """
    This will fetch an Ensembl specific primary id. This is used for extenral
    id in our database. It will be the transcript id with the version stripped,
    or standard name.
    """

    primary = transcript(feature)
    if not primary:
        primary = embl.standard_name(feature)
    assert primary, "Could not generate primary id for %s" % feature
    return primary.split('.', 1)[0]


def accession(feature):
    """
    This will compute the accession for the given feature. This will either be
    the transcript id or the standard name.
    """

    acc = transcript(feature) or embl.standard_name(feature)
    if not acc:
        raise ValueError("No accession possible for %s" % feature)
    return acc


def url(feature):
    """
    Get the URL for the given Enesmbl feature.
    """
    return URL.format(transcript=transcript(feature))


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


def note_data(feature):
    """
    This will parse the notes data of the feature to produce a dict of key
    value mappings. This will strip out the RNA type in the final dictonary.
    """

    notes_data = notes(feature)
    if len(notes_data) > 1:
        notes_data = notes_data[1:]
    grouped = embl.grouped_annotations(notes_data, '=')
    if 'transcript_id' not in grouped:
        grouped['transcript_id'] = [transcript(feature)]
    return grouped


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

    if raw_rna_type(feature) == 'scaRNA':
        return 'scaRNA'
    return None


def generate_related(entries):
    """
    This goes through all given entries, which are assumed to all be from the
    same gene, and thus splicing variants, and populates the related_sequences
    feature with the required related sequence information.
    """
    return data.related_isoforms(entries)
