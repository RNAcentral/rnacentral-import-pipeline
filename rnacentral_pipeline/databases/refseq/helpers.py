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

import attr

from rnacentral_pipeline.databases import data as dat

import rnacentral_pipeline.databases.helpers.embl as embl

URL = 'https://www.ncbi.nlm.nih.gov/nuccore/{primary_id}.{version}'

NCRNA = {
    'ncRNA',
    'precursor_RNA',
    'rRNA',
}


def optional_id(feature):
    for ref in feature.qualifiers['db_xref']:
        if ref.startswith('GeneID:'):
            return ref
    return None


def parent_accession(record):
    return record.annotations['accessions'][0]


def url(record):
    return URL.format(
        primary_id=primary_id(record),
        version=embl.seq_version(record),
    )


def description(record, feature):
    rna_type = embl.rna_type(feature).replace('_', ' ')
    name = record.description
    name = re.sub(r'\.$', '', name)
    name = re.sub(',.+$', ', ' + rna_type, name)
    return name


def primary_id(record):
    return record.annotations['accessions'][0]


def accession(record, feature):
    return '{primary_id}.{version}:{start}..{stop}:{feature_type}'.format(
        primary_id=primary_id(record),
        version=embl.seq_version(record),
        start=feature.location.start + 1,
        stop=feature.location.end,
        feature_type=feature.type,
    )


def xref_data(feature):
    feature_annotations = feature.qualifiers.get('db_xref', [])
    return embl.grouped_annotations(feature_annotations, ':')


def as_entry(record, source, feature):
    """
    Create an Entry based upon the record, source feature and ncRNA feature.
    """

    try:
        optional = optional_id(feature)
    except:
        optional = optional_id(source)

    return dat.Entry(
        primary_id=primary_id(record),
        accession=accession(record, feature),
        ncbi_tax_id=embl.taxid(record),
        database='REFSEQ',
        sequence=embl.sequence(record, feature),
        regions=[],
        rna_type=embl.rna_type(feature),
        url=url(record),
        seq_version=embl.seq_version(record),
        optional_id=optional,
        note_data={},
        xref_data=xref_data(feature),
        species=embl.species(record),
        common_name=embl.common_name(record),
        lineage=embl.lineage(record),
        gene=embl.gene(feature),
        locus_tag=embl.locus_tag(feature),
        product=embl.product(feature),
        parent_accession=parent_accession(record),
        project=embl.project(record),
        keywords=embl.keywords(record),
        organelle=embl.organelle(source),
        experiment=embl.experiment(feature),
        inference=embl.inference(feature),
        old_locus_tag=embl.old_locus_tag(feature),
        operon=embl.operon(feature),
        standard_name=embl.standard_name(feature),
        description=description(record, feature),
        mol_type=embl.mol_type(source),
        is_composite='N',
        gene_synonyms=embl.gene_synonyms(feature),
        references=embl.references(record),
    )


def ncrna_features(features):
    return [f for f in features if f.type in NCRNA]


def generate_related(entries):
    """
    This goes through all given entries, which are assumed to all be from the
    same gene, and thus splicing variants, and populates the related_sequences
    feature with the required related sequence information.
    """

    for first in entries:
        related = first.related_sequences
        for second in entries:
            if first == second:
                continue

            relationship = 'isoform'
            if first.rna_type == 'precursor_RNA':
                if second.rna_type == 'miRNA':
                    relationship = 'mature_product'
                else:
                    raise ValueError("Unknown type of relationship")

            elif first.rna_type == 'miRNA':
                if second.rna_type == 'precursor_RNA':
                    relationship = "precursor"
                elif second.rna_type == 'miRNA':
                    continue
                else:
                    raise ValueError("Unknown tyoe of relationship")

            related.append(dat.RelatedSequence(
                sequence_id=second.accession,
                relationship=relationship,
                coordinates=[],
                evidence=dat.RelatedEvidence.empty()
            ))
        yield attr.evolve(first, related_sequences=related)
