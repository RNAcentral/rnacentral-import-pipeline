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

import re
import attr

import luigi

from ensembl import data
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
        return 'OTTT' in entry.xref_data

    def method_for(self, instance, field):
        name = field.name
        if instance.database == 'GENCODE' and \
                hasattr(instance, 'gencode_' + field.name):
            name = 'gencode_' + field.name
        if hasattr(self, name):
            return getattr(self, name)
        return None

    def initial_entries(self, *args):
        parent = super(Gencode, self).initial_entries(*args)
        if self.is_gencode(parent[0]):
            parent.append(attr.assoc(
                parent[0],
                database='GENCODE',
            ))
        return parent

    def gencode_xref_data(self, summary, feature, current):
        xrefs = {k: v for k, v in current.xref_data.items() if k != 'OTTT'}
        xrefs['Ensembl'] = current.accession
        return xrefs

    def gencode_accession(self, summary, feature, current):
        return current.primary_id

    def gencode_primary_id(self, summary, feature, current):
        for value in current['db_xrefs']['OTTT']:
            if re.match(r'^OTT\w+T\d+$', value):
                return value
        raise ValueError("Cannot find GENCODE primary id")

    def gencode_description(self, summary, feature, current):
        name = current.species
        rna_type = current.ncrna_class
        gencode_id = current.primary_id
        assert name, "Must have a name to create description"
        assert rna_type, "Must have an rna_type to create description"
        return '{species} {rna_type} {gencode_id}'.format(
            species=name,
            rna_type=rna_type,
            gencode_id=gencode_id,
        )

    def gencode_references(self, summary, feature, current):
        return [data.Reference(
            authors=(
                "Jennifer Harrow, Adam Frankish, Jose M. Gonzalez, "
                "Electra Tapanari, Mark Diekhans, Felix Kokocinski, "
                "Bronwen L. Aken, Daniel Barrell, Amonida Zadissa, "
                "Stephen Searle, If Barnes, Alexandra Bignell, "
                "Veronika Boychenko, Toby Hunt, Mike Kay, Gaurab Mukherjee, "
                "Jeena Rajan, Gloria Despacio-Reyes, Gary Saunders, "
                "Charles Steward, Rachel Harte, Michael Lin, Cedric Howald, "
                "Andrea Tanzer, Thomas Derrien, Jacqueline Chrast, "
                "Nathalie Walters, Suganthi Balasubramanian, Baikang Pei, "
                "Michael Tress, Jose Manuel Rodriguez, Iakes Ezkurdia, "
                "Jeltje van Baren, Michael Brent, David Haussler, "
                "Manolis Kellis, Alfonso Valencia, Alexandre Reymond, "
                "Mark Gerstein, Roderic Guigo and Tim J. Hubbard"
            ),
            location="Genome Research",
            title=(
                "GENCODE: the reference human genome annotation for "
                "The ENCODE Project"
            ),
            pmid=22955987,
            doi="10.1101/gr.135350.111",
            accession=current.accession,
        )]

    def gencode_optional_id(self, *args):
        return None

    # def gencode_seq_version(self, *args):
    #     return ''


if __name__ == '__main__':
    luigi.run(main_task_cls=Gencode)
