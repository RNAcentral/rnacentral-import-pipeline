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

__doc__ = """
This is the module containing the importer for GENCODE data. GENCODE data is a
bit special in that it also has GENCODE annotations. These have to be handled
differently. Notably, we provide a separate xref for these annotations as well.
"""

import copy

import luigi

from ensembl.generic import EnsemblImporter as BaseImporter


class Gencode(BaseImporter):
    """
    Importer for GENCODE data from Ensembl. This handles the issues that come
    from having GENCODE as well as Ensembl data in the same file. Basically
    this will create 2 entries for any feature that comes from GENCODE. The
    first will be the standard Ensembl one, while the second will be to
    indicate the data comes from GENCODE as well.
    """

    def is_gencode(self, entry):
        """
        Check if this entry is also from gencode. We need to separate out these
        two databases despite being in the same database to provide users clear
        information about the high quality annotations.
        """
        return any(xref.startswith('OTTT:') for xref in entry['db_xrefs'])

    def gencode_references(self, _):
        return [{
            'authors': "Jennifer Harrow, Adam Frankish, Jose M. Gonzalez, Electra Tapanari, Mark Diekhans, Felix Kokocinski, Bronwen L. Aken, Daniel Barrell, Amonida Zadissa, Stephen Searle, If Barnes, Alexandra Bignell, Veronika Boychenko, Toby Hunt, Mike Kay, Gaurab Mukherjee, Jeena Rajan, Gloria Despacio-Reyes, Gary Saunders, Charles Steward, Rachel Harte, Michael Lin, Cédric Howald, Andrea Tanzer, Thomas Derrien, Jacqueline Chrast, Nathalie Walters, Suganthi Balasubramanian, Baikang Pei, Michael Tress, Jose Manuel Rodriguez, Iakes Ezkurdia, Jeltje van Baren, Michael Brent, David Haussler, Manolis Kellis, Alfonso Valencia, Alexandre Reymond, Mark Gerstein, Roderic Guigóand Tim J. Hubbard",
            'location': "Genome Research",
            'title': "GENCODE: the reference human genome annotation for The ENCODE Project",
            'pmid': 22955987,
            'doi': "10.1101/gr.135350.111",
        }]

    def gencode_xrefs(self, entry):
        xrefs = [v for v in entry['db_xrefs'] if not v.startswith('OTTT:')]
        xrefs.append('Ensembl:%s' % entry['accession'])
        return xrefs

    def gencode_accession(self, entry):
        return '{parent}:{gencode}:{type}'.format(
            parent=entry['parent_accession'],
            gencode=self.gencode_primary_id(entry),
            type=entry['ncrna_class'],
        )

    def gencode_primary_id(self, entry):
        for xref in entry['db_xrefs']:
            if xref.startswith('OTTT:OTTHUMT'):
                return xref[4:]
        raise ValueError("Cannot find GENCODE primary id")

    def gencode_entry(self, entry):
        """
        Create a new RNAcentralEntry for the GENCODE data in this feature. This
        assumes that the entry has GENCODE data.
        """

        return entry.update({
            'primary_id': self.gencode_primary_id(entry),
            'accession': self.gencode_accession(entry),
            'db_xrefs': self.gencode_xrefs(entry),
            'references': self.gencode_references(entry),
        })

    def rnacentral_entries(self, annotations, feature,
                           ignore_nongencode=False, gencode_only=False,
                           **kwargs):
        """
        Compute the RNAcentralEntrys for the feature. The features
        """

        entries = super(Gencode, self).rnacentral_entries(annotations, feature)
        for entry in entries:
            if self.is_gencode(entry):
                if not gencode_only:
                    yield entry
                yield self.gencode_entry(copy.deepcopy(entry))
            elif not ignore_nongencode and not gencode_only:
                yield entry


if __name__ == '__main__':
    luigi.run(main_task_cls=Gencode)
