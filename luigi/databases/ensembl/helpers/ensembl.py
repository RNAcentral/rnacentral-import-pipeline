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

from databases.data import Reference

from ..helpers import bio as helpers
from ..rna_type_inference import RnaTypeInference

URL = 'http://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript}'


def exons(feature):
    """
    Determine all Exons in this feature.
    """

    parts = [feature.location]
    if hasattr(feature.location, 'parts'):
        parts = feature.location.parts
    return [helpers.exon(l) for l in parts]


def references(accession):
    """
    Get the standard reference for all Ensembl entries.
    """
    return [Reference(
        accession=accession,
        authors=(
            "Aken BL, Ayling S, Barrell D, Clarke L, Curwen V, Fairley "
            "S, Fernandez Banet J, Billis K, Garci a Giro n C, Hourlier "
            "T, Howe K, Kahari A, Kokocinski F, Martin FJ, Murphy DN, "
            "Nag R, Ruffier M, Schuster M, Tang YA, Vogel JH, White "
            "S, Zadissa A, Flicek P, Searle SM."
        ),
        location="Database (Oxford). 2016 Jun 23",
        title="The Ensembl gene annotation system",
        pmid=27337980,
        doi="10.1093/database/baw093",
    )]


def is_pseudogene(summary, feature):
    """
    Check if the given feature is a pseudogene. This is determined by check
    if the feature is labeled as a psuedogene or the gene it is a transcript
    from a pseudogene. """

    notes = feature.qualifiers.get('note', [])
    if any('pseudogene' in n for n in notes):
        return True

    gene = helpers.gene(feature)
    return summary.is_pseudogene(gene)


def rna_type(feature, xref_data):
    """
    Compute the RNA type of the given feature.
    """
    base_type = helpers.rna_type(feature)
    inference = RnaTypeInference()
    return inference.infer_rna_type(xref_data, base_type)


def feature_description(summary, feature, entry):
    """
    Create a description using the feature for the given entry.
    """
    trimmed = summary.trimmed_description(feature)
    if not trimmed:
        return None

    return '{species} {trimmed}'.format(
        species=entry.species,
        trimmed=trimmed,
    )


def locus_description(entry):
    """
    Generate a description for the entry based upon the locus this is a part
    of. This will be of the form 'Homo sapiens (human) lncRNA xist
    """
    species = entry.species
    if entry.common_name:
        species += ' (%s)' % entry.common_name

    assert entry.rna_type, "Cannot build description without rna_type"
    return '{species} {rna_type} {locus_tag}'.format(
        species=species,
        rna_type=entry.rna_type.replace('_', ' '),
        locus_tag=entry.locus_tag or '',
    ).strip()


def description(summary, feature, entry):
    """
    Generate a description of rthe given entry. This will use either the
    feature description or the locus description if that required.
    """
    return feature_description(summary, feature, entry) or \
        locus_description(entry)


def seq_version(feature):
    """
    Compute the sequence version, if any, of the given id.
    """

    acc = accession(feature)
    if '.' in acc:
        return acc.split('.', 1)[1]
    return None


def primary_id(feature):
    """
    This will fetch an Ensembl specific primary id. This is used for extenral
    id in our database. It will be the transcript id with the version stripped,
    or standard name.
    """

    primary = helpers.transcript(feature)
    if not primary:
        primary = helpers.standard_name(feature)
    assert primary, "Could not generate primary id for %s" % feature
    return primary.split('.', 1)[0]


def accession(feature):
    """
    This will compute the accession for the given feature. This will either be
    the transcript id or the standard name.
    """

    acc = helpers.transcript(feature) or helpers.standard_name(feature)
    if not acc:
        raise ValueError("No accession possible for %s" % feature)
    return acc


def url(feature):
    """
    Get the URL for the given Enesmbl feature.
    """
    return URL.format(transcript=helpers.transcript(feature))
