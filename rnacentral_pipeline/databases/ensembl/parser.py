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

import attr

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import embl

from . import helpers
from .data import Context


def as_entry(record, gene, feature, context):
    """
    Turn the Record, Gene feature, transcript feature and Context into a Entry
    object for output.
    """

    species, common_name = helpers.organism_naming(record)
    sequence = str(feature.extract(record.seq))

    xref_data = embl.xref_data(feature)
    accession = helpers.accession(feature)
    chromosome = helpers.chromosome(record)
    exons = helpers.exons(record, feature)

    entry = data.Entry(
        primary_id=helpers.primary_id(feature),
        accession=accession,
        ncbi_tax_id=embl.taxid(record),
        database='ENSEMBL',
        sequence=sequence,
        exons=exons,
        rna_type=helpers.rna_type(context.inference, feature, xref_data),
        url=helpers.url(feature),
        seq_version=helpers.seq_version(feature),
        lineage=embl.lineage(record),
        chromosome=chromosome,
        parent_accession=record.id,
        common_name=common_name,
        species=species,
        gene=embl.locus_tag(gene),
        locus_tag=embl.locus_tag(gene),
        optional_id=gene,
        note_data=helpers.note_data(feature),
        xref_data=xref_data,
        product=helpers.product(feature),
        references=helpers.references(),
        mol_type='genomic DNA',
        pseudogene='N',
        is_composite='N',
    )

    return attr.assoc(
        entry,
        description=helpers.description(context, feature, entry)
    )


def parse(raw, family_file):
    """
    This will parse an EMBL file for all Ensembl Entries to import.
    """

    context = Context(family_file)
    for (record, gene, feature) in embl.transcripts(raw):
        if context.is_supressed(feature):
            continue
        yield as_entry(record, gene, feature, context)
