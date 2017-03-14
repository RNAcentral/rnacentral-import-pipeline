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

import luigi

from rnacentral_entry import RNAcentralEntry

from ensembl.base import BioImporter
from ensembl.base import qualifier_value


class EnsemblImporter(BioImporter):
    """
    The generic Ensembl data importer. This mostly does the right thing for all
    Ensembl data, however, there are some specific cases where it does not.
    Notably, some organisms like Mouse and Human have GENCODE data as well. In
    those cases a more specific importer must be used to get the correct data.
    """

    input_file = luigi.Parameter()
    test = luigi.BoolParameter(default=False, significant=False)
    destination = luigi.Parameter(default='/tmp')

    def format(self):
        """
        Determine the format. This is always 'embl' for Ensembl.
        """
        return 'embl'

    def output(self):
        """
        Create all luigi output objects.
        """
        return self.build_output(self.destination, self.input_file)

    def assembly_info(self, feature):
        def as_exon(location):
            """
            Turn a biopython location into the dict we use.
            """

            return {
                'primary_start': location.start + 1,
                'primary_end': int(location.end),
                'complement': location.strand == -1,
            }

        parts = [feature.location]
        if hasattr(feature.location, 'parts'):
            parts = feature.location.parts
        return [as_exon(l) for l in parts]

    def references(self, _):
        return [{
            'authors': "Andrew Yates, Wasiu Akanni, M. Ridwan Amode, Daniel Barrell, Konstantinos Billis, Denise Carvalho-Silva, Carla Cummins, Peter Clapham, Stephen Fitzgerald, Laurent Gil1 Carlos Garcín Girón, Leo Gordon, Thibaut Hourlier, Sarah E. Hunt, Sophie H. Janacek, Nathan Johnson, Thomas Juettemann, Stephen Keenan, Ilias Lavidas, Fergal J. Martin, Thomas Maurel, William McLaren, Daniel N. Murphy, Rishi Nag, Michael Nuhn, Anne Parker, Mateus Patricio, Miguel Pignatelli, Matthew Rahtz, Harpreet Singh Riat, Daniel Sheppard, Kieron Taylor, Anja Thormann, Alessandro Vullo, Steven P. Wilder, Amonida Zadissa, Ewan Birney, Jennifer Harrow, Matthieu Muffato, Emily Perry, Magali Ruffier, Giulietta Spudich, Stephen J. Trevanion, Fiona Cunningham, Bronwen L. Aken, Daniel R. Zerbino, Paul Flicek",
            'location': "Nucleic Acids Res. 2016 44 Database issue:D710-6",
            'title': "Ensembl 2016",
            'pmid': 26687719,
            'doi': "10.1093/nar/gkv115",
        }]

    def standard_annotations(self, record):
        pattern = re.compile(r'\((.+)\)$')
        common_name = None
        species = record.annotations['organism']
        match = re.search(pattern, species)
        if match:
            common_name = match.group(1)
            species = re.sub(pattern, '', species).strip()

        taxid = None
        source = record.features[0]
        if source.type == 'source':
            taxid = int(qualifier_value(source, 'db_xref', r'^taxon:(\d+)$'))

        return {
            'accession': None,
            'database': 'ENSEMBL',
            'lineage': '; '.join(record.annotations['taxonomy']),
            'mol_type': 'genomic DNA',
            'pseudogene': 'N',
            'parent_accession': record.id,
            'seq_version': '',
            'common_name': common_name,
            'species': species,
            'ncbi_tax_id': taxid,
            'is_composite': 'N',
            'references': self.references(record),
            'sequence': record.seq,
        }

    def description(self, annotations, feature):
        species = annotations['species']
        if annotations.get('common_name', None):
            species += ' (%s)' % annotations['common_name']

        rna_type = self.ncrna(feature) or feature.type
        return '{species} {rna_type}'.format(
            species=species,
            rna_type=rna_type,
        )

    def gene(self, feature):
        return qualifier_value(feature, 'gene', '^(.+)$')

    def transcript(self, feature):
        return qualifier_value(feature, 'note', '^transcript_id=(.+)$')

    def ncrna(self, feature):
        return feature.qualifiers.get('note', [''])[0]

    def note(self, feature):
        note = feature.qualifiers.get('note', [])
        if len(note) > 1:
            return note[1:]
        return None

    def primary_id(self, annotations, feature):
        transcript = self.transcript(feature)
        ncrna = self.ncrna(feature)
        assert transcript, 'Cannot create a primary id without transcript id'
        assert ncrna, 'Cannot create a primary id without ncRNA type'
        assert annotations['parent_accession']

        return '{parent}:{transcript}:{type}'.format(
            parent=annotations['parent_accession'],
            transcript=transcript,
            type=ncrna,
        )

    def entry_specific_data(self, feature):
        start, end = sorted([int(feature.location.start),
                             int(feature.location.end)])

        # Ensure that this is 1 based, even though biopython
        # switches to 0 based
        start += 1

        return {
            'primary_id': self.transcript(feature),
            'assembly_info': self.assembly_info(feature),
            'db_xrefs': feature.qualifiers.get('db_xref', []),
            'feature_location_start': start,
            'feature_location_end': end,
            'feature_type': feature.type,
            'gene': self.gene(feature),
            'ncrna_class': self.ncrna(feature),
            'note': self.note(feature),
            'locus_tag': feature.qualifiers.get('locus_tag', ''),
        }

    def rnacentral_entries(self, annotations, feature, **kwargs):
        data = dict(annotations)
        data.update(self.entry_specific_data(feature))
        data['accession'] = self.primary_id(annotations, feature)
        data['description'] = self.description(annotations, feature)
        data['sequence'] = self.sequence(annotations, feature)
        return [data]


if __name__ == '__main__':
    luigi.run(main_task_cls=EnsemblImporter)
