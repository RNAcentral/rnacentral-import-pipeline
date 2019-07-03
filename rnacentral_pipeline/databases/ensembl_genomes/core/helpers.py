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

import attr

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.helpers import publications as pubs

from rnacentral_pipeline.databases.ensembl import helpers as ensembl

EXCLUDED_TYPES = {
    'nontranslating_CDS',
    'sense_intronic',
    'transposable_element',
}

NORMALIZED_RNA_TYPES = {
    'pre_miRNA': 'precursor_RNA',
    'atlncRNA': 'lncRNA',
    'sense_intronic': 'lncRNA',
    'atRNA': 'antisense_RNA',
    'antisense': 'antisense_RNA',
    'otherRNA': 'other',
    'lincRNA': 'lncRNA',
}


def is_pseudogene(gene, feature):
    raw_notes = feature.qualifiers.get('note', [])
    if gene:
        raw_notes.extend(gene.qualifiers.get('note', []))
    return any('pseudogene' in n for n in raw_notes)


def is_ncrna(feature):
    if not ensembl.is_ncrna(feature):
        return False
    rna_type = ensembl.raw_rna_type(feature)
    return rna_type not in EXCLUDED_TYPES


def primary_id(feature):
    """
    This will fetch an Ensembl specific primary id. This is used for extenral
    id in our database. It will be the transcript id with the version stripped,
    or standard name.
    """

    primary = ensembl.transcript(feature)
    if not primary:
        primary = embl.standard_name(feature)
    assert primary, "Could not generate primary id for %s" % feature
    return primary


def rna_type(feature):
    raw = ensembl.raw_rna_type(feature)
    return NORMALIZED_RNA_TYPES.get(raw, raw)


def seq_version(feature):
    return ensembl.seq_version(feature) or '1'


def description(gene, entry):
    species = entry.species
    if entry.common_name:
        species += ' (%s)' % entry.common_name

    gene_name = ensembl.notes(gene)
    if gene_name and len(gene_name) == 1:
        gene_name = gene_name[0]
        if gene_name.endswith(']'):
            gene_name = re.sub(r'\s*\[.+\]$', '', gene_name)

        gene_name.strip()
        if not re.search('^U\d+', gene_name):
            gene_name = gene_name[0].lower() + gene_name[1:]

        return '{species} {gene_name}'.format(
            species=species,
            gene_name=gene_name
        )

    assert entry.rna_type, "Cannot build description without rna_type"
    return '{species} {rna_type} {locus_tag}'.format(
        species=species,
        rna_type=entry.rna_type.replace('_', ' '),
        locus_tag=entry.locus_tag or entry.gene or '',
    ).strip()


def xref_data(feature):
    result = {}
    for key, values in embl.xref_data(feature).items():
        if key == 'RefSeq_dna':
            result['RefSeq'] = values
        elif key == 'TAIR_LOCUS_MODEL':
            result['TAIR'] = values
        else:
            result[key] = values
    return result


def as_entry(context, record, current_gene, feature):
    species, common_name = ensembl.organism_naming(record)

    entry = data.Entry(
        primary_id=primary_id(feature),
        accession=context.accession(primary_id(feature)),
        ncbi_tax_id=embl.taxid(record),
        database=context.database,
        sequence=embl.sequence(record, feature),
        regions=ensembl.regions(record, feature),
        rna_type=rna_type(feature),
        url='',
        seq_version=seq_version(feature),
        lineage=embl.lineage(record),
        common_name=common_name,
        species=species,
        gene=embl.gene(current_gene),
        locus_tag=embl.locus_tag(current_gene),
        xref_data=xref_data(feature),
        product=ensembl.product(feature),
        references=context.references,
    )

    return attr.evolve(
        entry,
        description=description(current_gene, entry)
    )
