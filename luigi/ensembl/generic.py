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

import luigi

from ensembl.base import BioImporter
from ensembl import data
from ensembl import helpers
from ensembl.rna_type_inference import RnaTypeInference
from rfam import utils as rfutils

import parameters

LOGGER = logging.getLogger(__name__)


class EnsemblImporter(BioImporter):
    """
    The generic Ensembl data importer. This mostly does the right thing for all
    Ensembl data, however, there are some specific cases where it does not.
    Notably, some organisms like Mouse and Human have GENCODE data as well. In
    those cases a more specific importer must be used to get the correct data.
    """

    input_file = parameters.GenericFileParameter()
    test = luigi.BoolParameter(default=False, significant=False)
    destination = parameters.PathParameter(default='/tmp')

    def format(self):
        """
        Determine the format. This is always 'embl' for Ensembl.
        """
        return 'embl'

    def is_from_suppressed_rfam_model(self, current):
        """
        Check fi the currenty entry is from a suppressed Rfam model. Some Rfam
        models are not good fits for the RNAcentral, notably the models of part
        of an lncRNA or piRNA models. We can compute which models are not good
        fits and so we use that to reject some Ensembl annotation as they are
        the result of searching using one of these models.

        :param data.Entry entry: The Entry to check.
        :return bool: True if this is entry is from a suppressed Rfam model.
        """
        inference = RnaTypeInference()
        rfam_model = inference.rfam_xref(current)
        if not rfam_model:
            return False
        name = inference.rfam_name(rfam_model)
        if name is None:
            return False

        mapping = rfutils.name_to_suppression()
        if name not in mapping:
            raise ValueError("Unknown Rfam model name: %s" % name)
        return mapping[name]

    def initial_entries(self, record, summary, feature):
        try:
            sequence = str(feature.extract(record.seq))
        except Exception as err:
            LOGGER.exception(err)
            LOGGER.warn("Could not get sequence for %s", feature)
            return []

        species, common_name = helpers.organism_naming(record)
        transcript_id = helpers.transcript(feature)
        gene = helpers.gene(feature)

        primary_id = transcript_id
        accession = transcript_id
        standard_name = helpers.standard_name(feature)

        if not transcript_id:
            primary_id = standard_name
            accession = standard_name

        entry = data.Entry(
            primary_id=primary_id,
            accession=accession,
            seq=sequence,
            ncbi_tax_id=helpers.taxid(record),
            database='ENSEMBL',
            lineage=helpers.lineage(record),
            chromosome=helpers.chromosome(record),
            parent_accession=record.id,
            common_name=common_name,
            species=species,
            gene=gene,
            locus_tag=summary.locus_tag(gene),
            optional_id=gene,
            note_data=helpers.note_data(feature),
            xref_data=helpers.xref_data(feature),
        )

        if self.is_from_suppressed_rfam_model(entry):
            LOGGER.debug("Skipping feature %s because it is from a suppressed"
                         " Rfam family", feature)
            return []

        return [entry]

    def description(self, summary, feature, current):
        super_method = super(EnsemblImporter, self).description
        computed = super_method(summary, feature, current)
        if computed:
            return computed

        species = current.species
        if current.common_name:
            species += ' (%s)' % current.common_name

        assert current.rna_type, "Cannot build description without rna_type"
        return '{species} {rna_type} {locus_tag}'.format(
            species=species,
            rna_type=current.rna_type,
            locus_tag=current.locus_tag,
        )

    def exons(self, summary, feature, current):
        parts = [feature.location]
        if hasattr(feature.location, 'parts'):
            parts = feature.location.parts
        return [data.Exon.from_biopython(l) for l in parts]

    def references(self, summary, feature, current):
        return [data.Reference(
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
            accession=current.accession,
        )]

    def seq_version(self, current, *args):
        transcript = current.primary_id
        if '.' in transcript:
            parts = transcript.split('.', 1)
            return parts[1]
        return ''

    def product(self, summary, feature, current):
        if helpers.rna_type(feature) == 'scaRNA':
            return 'scaRNA'
        return ''

    def rna_type(self, summary, feature, current):
        inference = RnaTypeInference()
        base_type = helpers.rna_type(feature)
        found = inference.infer_rna_type(current, base_type)
        return found


if __name__ == '__main__':
    luigi.run(main_task_cls=EnsemblImporter)
