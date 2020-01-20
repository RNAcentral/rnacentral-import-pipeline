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

import logging

from Bio import SeqIO

from rnacentral_pipeline.databases.data import Entry

import rnacentral_pipeline.databases.helpers.embl as embl

from . import dr
from . import helpers
from . import mapping as tpa

LOGGER = logging.getLogger(__name__)


class InvalidEnaFile(Exception):
    """
    This is raised when there is something wrong with the ENA EMBL file.
    """
    pass


def parse(handle):
    """
    Parse a file like object into an iterable of Entry objects. This will parse
    each feature in all records of the given EMBL formatted file to produce the
    Entry objects.
    """

    dr_mapping = dr.mapping(handle)
    handle.seek(0)
    for record in SeqIO.parse(handle, 'embl'):

        if len(record.features) != 2:
            raise InvalidEnaFile("ENA EMBL files must have 2 features/record")

        if helpers.is_protein(record.features[1]):
            LOGGER.info("Skipping mis-annotated protein: %s", record.id)
            continue

        if helpers.is_pseudogene(record.features[1]):
            LOGGER.info("Skipping pseudogene")
            continue

        if record.id not in dr_mapping:
            raise InvalidEnaFile("Somehow parsed DR refs are for wrong record")

        record_refs = dr_mapping[record.id]
        accession = helpers.accession(record)

        feature = record.features[1]
        prod = helpers.product(feature)
        if prod:
            prod = prod[0:500]

        entry = Entry(
            primary_id=helpers.primary_id(feature),
            accession=accession,
            ncbi_tax_id=helpers.taxid(record),
            database='ENA',
            sequence=helpers.sequence(record),
            regions=[],
            rna_type=helpers.rna_type(feature),
            url=helpers.url(record),
            seq_version=embl.seq_version(record),

            note_data=helpers.note_data(feature),
            xref_data=helpers.xref_data(record, feature, record_refs),

            chromosome=helpers.chromosome(record),

            species=helpers.species(record),
            common_name=helpers.common_name(record),
            lineage=helpers.lineage(record),

            gene=embl.gene(feature),
            locus_tag=embl.locus_tag(feature),
            product=prod,
            parent_accession=helpers.parent_accession(record),
            project=embl.project(record),
            keywords=helpers.keywords(record),
            organelle=helpers.organelle(record),
            anticodon=helpers.anticodon(record, feature),
            experiment=embl.experiment(feature),
            function=helpers.function(feature),
            inference=embl.inference(feature),
            old_locus_tag=embl.old_locus_tag(feature),
            operon=helpers.operon(feature),
            standard_name=embl.standard_name(feature),
            description=helpers.description(record),
            mol_type=helpers.mol_type(record),
            is_composite=helpers.is_composite(feature),

            gene_synonyms=helpers.gene_synonyms(feature),
            references=helpers.references(record, feature),
        )

        if helpers.is_skippable_sequence(entry):
            continue

        yield entry


def parse_with_mapping_file(handle, mapping_handle):
    mapping = tpa.load(mapping_handle)
    mapping.validate()
    return tpa.apply(mapping, parse(handle))
