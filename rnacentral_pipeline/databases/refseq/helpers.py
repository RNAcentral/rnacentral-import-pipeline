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

from rnacentral_pipeline.databases.data import Entry

import rnacentral_pipeline.databases.helpers.embl as embl

URL = 'https://www.ncbi.nlm.nih.gov/nuccore/{primary_id}.{version}'


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


def description(record):
    return record.description


def primary_id(record):
    return record.annotations['accessions'][0]


def accession(record, feature):
    return '{primary_id}.{version}:{start}..{stop}:{feature_type}'.format(
        versioned=primary_id(record),
        version=embl.seq_version(record),
        start=feature.start,
        stop=feature.stop,
        feature_type=feature.type,
    )


def as_entry(record, source, feature):
    """
    Modify an ENA entry into an approbate RefSeq entry.
    """

    # return attr.evolve(
    #     parent_accession=parent_accession(entry),
    # )

    return Entry(
        primary_id=primary_id(record),
        accession=accession(record, feature),
        ncbi_tax_id=embl.taxid(record),
        database='REFSEQ',
        sequence=sequence(record),
        exons=[],
        rna_type=embl.rna_type(feature),
        url=url(record),
        seq_version=embl.seq_version(record),

        optional_id=optional_id(record),

        note_data=ena.note_data(feature),
        xref_data=xref_data(record, feature, record_refs),

        chromosome=embl.chromosome(source),

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
        description=description(record),
        mol_type=embl.mol_type(source),
        is_composite='N',

        gene_synonyms=embl.gene_synonyms(feature),
        references=embl.references(record),
    )
