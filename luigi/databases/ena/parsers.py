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


def parse(handle):
    """
    Parse a file like object into an iterable of Entry objects. This will parse
    each feature in all records of the given EMBL formatted file to produce the
    Entry objects.
    """

    for record in SeqIO.parse(handle, 'embl'):
        taxid = embl.taxid(record)
        species = embl.species(record)
        chromosome = helpers.chromosome(record)
        common_name = embl.common_name(record)
        lineage = embl.lineage(record)
        project = embl.project(record)
        description = embl.description(record)
        seq_version = embl.seq_version(record)
        accession = helpers.accession(record)
        division = embl.division(record)
        references = embl.references(accession, record)
        sequence = helpers.sequence(record)

        for feature in record.features[1:]:
            yield Entry(
                primary_id=helpers.primary_id(feature),
                accession=accession,
                ncbi_tax_id=taxid,
                database='ENA',
                sequence=sequence,
                exons=helpers.exons(feature),
                rna_type=helpers.rna_type(feature),
                url=helpers.url(feature),
                seq_version=seq_version,
                note_data=helpers.note_data(feature),
                xref_data=embl.xref_data(feature),
                chromosome=chromosome,
                species=species,
                common_name=common_name,
                lineage=lineage,
                gene=embl.gene(feature),
                locus_tag=embl.locus_tag(feature),
                product=helpers.product(feature),
                parent_accession=helpers.parent_accession(record),
                ordinal=helpers.ordinal(feature),
                non_coding_id=helpers.non_coding_id(feature),
                project=project,
                keywords=helpers.keywords(record),
                division=division,
                organelle=helpers.organelle(feature),
                allele=helpers.allele(feature),
                anticodon=helpers.anticodon(feature),
                experiment=helpers.experiment(feature),
                function=helpers.function(feature),
                inference=helpers.inference(feature),
                map=helpers.map(feature),
                old_locus_tag=helpers.old_locus_tag(feature),
                operon=helpers.operon(feature),
                standard_name=embl.standard_name(feature),
                description=description,
                mol_type=helpers.mol_type(record),
                is_composite=helpers.is_composite(feature),
                pseudogene=helpers.pseudogene(feature),
                gene_synonyms=helpers.gene_synonyms(feature),
                references=references,
            )
