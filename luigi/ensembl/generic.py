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

import luigi

from ensembl.base import BioImporter
from ensembl import data
from ensembl import helpers
from ensembl.rna_type_inference import RnaTypeInference

import parameters


class EnsemblImporter(BioImporter):
    """
    The generic Ensembl data importer. This mostly does the right thing for all
    Ensembl data, however, there are some specific cases where it does not.
    Notably, some organisms like Mouse and Human have GENCODE data as well. In
    those cases a more specific importer must be used to get the correct data.
    """

    input_file = parameters.FileParameter()
    test = luigi.BoolParameter(default=False, significant=False)
    destination = luigi.Parameter(default='/tmp')

    def format(self):
        """
        Determine the format. This is always 'embl' for Ensembl.
        """
        return 'embl'

    def initial_entries(self, record, summary, feature):
        species, common_name = helpers.organism_naming(record)
        transcript_id = helpers.transcript(feature)
        gene = helpers.gene(feature)

        return [data.Entry(
            primary_id=transcript_id,
            accession=transcript_id,
            seq=str(feature.extract(record.seq)),
            ncbi_tax_id=helpers.taxid(record),
            database='ENSEMBL',
            lineage=helpers.lineage(record),
            chromosome=record.id.split('.')[0],
            parent_accession=record.id,
            common_name=common_name,
            species=species,
            gene=gene,
            optional_id=gene,
            note_data=helpers.note_data(feature),
            locus_tag=helpers.locus_tag(feature),
            xref_data=helpers.xref_data(feature),
        )]

    def description(self, summary, feature, current):
        if current.gene in summary.gene_info:
            super_method = super(EnsemblImporter, self).description
            computed = super_method(summary, feature, current)
            if computed:
                return computed

        species = current.species
        if current.common_name:
            species += ' (%s)' % current.common_name

        rna_type = current.rna_type
        transcript = current.primary_id
        assert rna_type, "Cannot build description without rna_type"
        assert transcript, "Cannot build description without transcript"

        return '{species} {rna_type} transcript {transcript}'.format(
            species=species,
            rna_type=rna_type,
            transcript=transcript,
        )

    def exons(self, summary, feature, current):
        parts = [feature.location]
        if hasattr(feature.location, 'parts'):
            parts = feature.location.parts
        return [data.Exon.from_biopython(l) for l in parts]

    def references(self, summary, feature, current):
        return [data.Reference(
            authors=(
                "Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, "
                "Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, "
                "Peter Clapham, Stephen Fitzgerald, Laurent Gil Carlos Garci"
                "n Giro n, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, "
                "Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, "
                "Stephen Keenan, Ilias Lavidas, Fergal J. Martin, "
                "Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, "
                "Michael Nuhn, Anne Parker, Mateus Patricio, "
                "Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, "
                "Daniel Sheppard, Kieron Taylor, Anja Thormann, "
                "Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, "
                "Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, "
                "Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, "
                "Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, "
                "Paul Flicek"
            ),
            location="Nucleic Acids Res. 2016 44 Database issue:D710-6",
            title="Ensembl 2016",
            pmid=26687719,
            doi="10.1093/nar/gkv115",
            accession=current.accession,
        )]

    def seq_version(self, current, *args):
        transcript = current.primary_id
        if '.' in transcript:
            parts = transcript.split('.', 1)
            return parts[1]
        return ''

    def product(self, summary, feature, current):
        if current.rna_type == 'scaRNA':
            return current.rna_type
        return ''

    def rna_type(self, summary, feature, current):
        inference = RnaTypeInference()
        return inference.infer_rna_type(current)


if __name__ == '__main__':
    luigi.run(main_task_cls=EnsemblImporter)
