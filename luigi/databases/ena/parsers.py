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

from Bio import SeqIO

from databases.data import Entry

from . import helpers
import databases.helpers.embl as embl


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

    for record in SeqIO.parse(handle, 'embl'):

        if len(record.features) != 2:
            raise InvalidEnaFile("ENA EMBL files must have 2 features/record")

        accession = helpers.accession(record)

        feature = record.features[1]
        yield Entry(
            primary_id=helpers.primary_id(feature),
            accession=accession,
            ncbi_tax_id=helpers.taxid(record),
            database='ENA',
            sequence=helpers.sequence(record),
            exons=helpers.exons(record, feature),
            rna_type=helpers.rna_type(feature),
            url=helpers.url(record),
            seq_version=embl.seq_version(record),

            note_data=helpers.note_data(feature),
            xref_data=embl.xref_data(feature),

            chromosome=helpers.chromosome(record),

            species=helpers.species(record),
            common_name=helpers.common_name(record),
            lineage=helpers.lineage(record),

            gene=embl.gene(feature),
            locus_tag=embl.locus_tag(feature),
            product=helpers.product(feature),
            parent_accession=helpers.parent_accession(record),
            ordinal=helpers.ordinal(feature),
            non_coding_id=helpers.non_coding_id(feature),
            project=embl.project(record),
            keywords=helpers.keywords(record),
            division=helpers.division(record),
            organelle=helpers.organelle(record),
            anticodon=helpers.anticodon(record, feature),
            experiment=embl.experiment(feature),
            function=helpers.function(record),
            inference=embl.inference(feature),
            old_locus_tag=helpers.old_locus_tag(feature),
            operon=helpers.operon(feature),
            standard_name=embl.standard_name(feature),
            description=embl.description(record),
            mol_type=helpers.mol_type(record),
            is_composite=helpers.is_composite(feature),
            pseudogene=helpers.pseudogene(feature),

            gene_synonyms=helpers.gene_synonyms(feature),
            references=embl.references(accession, record),
        )
