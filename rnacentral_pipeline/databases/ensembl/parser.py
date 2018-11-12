# -*- coding: utf-8 -*-

"""
Copyright [2010-2018] EMBL-European Bioinformatics Institute
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

import logging

import attr

from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.gencode import helpers as gencode

from . import helpers
from .data import Context

LOGGER = logging.getLogger(__name__)

IGNORE_FEATURES = {
    'source',
    'STS',
    'misc_feature',
}


def as_entry(record, gene, feature, context):
    """
    Turn the Record, Gene feature, transcript feature and Context into a Entry
    object for output.
    """

    species, common_name = helpers.organism_naming(record)
    xref_data = embl.xref_data(feature)

    try:
        sequence = embl.sequence(record, feature)
    except Exception as err:
        LOGGER.exception(err)
        return None

    entry = data.Entry(
        primary_id=helpers.primary_id(feature),
        accession=helpers.accession(feature),
        ncbi_tax_id=embl.taxid(record),
        database='ENSEMBL',
        sequence=sequence,
        regions=helpers.regions(record, feature),
        rna_type=helpers.rna_type(context.inference, feature, xref_data),
        url=helpers.url(feature),
        seq_version=helpers.seq_version(feature),
        lineage=embl.lineage(record),
        chromosome=helpers.chromosome(record),
        parent_accession=record.id,
        common_name=common_name,
        species=species,
        gene=embl.locus_tag(gene),
        locus_tag=embl.locus_tag(gene),
        optional_id=embl.gene(gene),
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
        description=helpers.description(gene, entry)
    )


def parse(raw, family_file, gencode_file=None):
    """
    This will parse an EMBL file for all Ensembl Entries to import.
    """

    context = Context.build(family_file, gencode_file=gencode_file)
    for record in SeqIO.parse(raw, 'embl'):
        current_gene = None
        ncrnas = []
        for feature in record.features:
            if feature.type in IGNORE_FEATURES:
                continue

            if embl.is_gene(feature):
                current_gene = feature
                for entry in helpers.generate_related(ncrnas):
                    yield entry
                    if context.from_gencode(entry):
                        yield gencode.update_entry(entry)
                ncrnas = []
                continue

            if helpers.is_pseudogene(current_gene, feature):
                continue

            if not helpers.is_ncrna(feature):
                continue

            entry = as_entry(record, current_gene, feature, context)
            if not entry or context.is_supressed(entry):
                continue

            ncrnas.append(entry)

        for entry in helpers.generate_related(ncrnas):
            yield entry
            if context.from_gencode(entry):
                yield gencode.update_entry(entry)
