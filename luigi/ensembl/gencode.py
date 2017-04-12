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
                hasattr(self, 'gencode_' + field.name):
            name = 'gencode_' + field.name
        if hasattr(self, name):
            return getattr(self, name)
        return None

    def initial_entries(self, *args):
        for entry in super(Gencode, self).initial_entries(*args):
            yield entry
            if self.is_gencode(entry):
                yield attr.assoc(entry, database='GENCODE')

    def gencode_xref_data(self, summary, feature, current):
        xrefs = {k: v for k, v in current.xref_data.items() if k != 'OTTT'}
        xrefs['Ensembl'] = [current.accession]
        return xrefs

    def gencode_primary_id(self, summary, feature, current):
        for value in current.xref_data['OTTT']:
            if re.match(r'^OTT\w+T\d+$', value):
                return value
        raise ValueError("Cannot find GENCODE primary id")

    def gencode_references(self, summary, feature, current):
        return [data.Reference(
            authors=(
                "Harrow J, Frankish A, Gonzalez JM, Tapanari E, Diekhans M, "
                "Kokocinski F, Aken BL, Barrell D, Zadissa A, Searle S, Barnes"
                " I, Bignell A, Boychenko V, Hunt T, Kay M, Mukherjee G, Rajan"
                " J,  Despacio-Reyes G, Saunders G, Steward C, Harte R, Lin M,"
                " Howald  C, Tanzer A, Derrien T, Chrast J, Walters N, "
                "Balasubramanian S,  Pei B, Tress M, Rodriguez JM, Ezkurdia "
                "I, van Baren J, Brent M,  Haussler D, Kellis M, Valencia A, "
                "Reymond A, Gerstein M, Guigo  R, Hubbard TJ. "
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

    def gencode_accession(self, summary, feature, current):
        return current.primary_id

    def gencode_optional_id(self, *args):
        return None


if __name__ == '__main__':
    luigi.run(main_task_cls=Gencode)
