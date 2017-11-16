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

from databases.ensembl.parsers import GencodeParser as Gencode

from tests.ensembl.utils import Base


class GencodeTests(Base):  # pylint: disable=R0904,C0111
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = Gencode

    def test_knows_if_entry_is_gencode(self):
        val = self.entries_from('GENCODE', "ENST00000535849.1")
        assert len(val) == 1

    def test_it_uses_correct_primary_id(self):
        entry = self.entry_from('GENCODE', "ENST00000400706.3")
        assert entry.primary_id == 'ENST00000400706'

    def test_it_keeps_xref_to_OTT(self):
        entry = self.entry_from('GENCODE', "ENST00000540226.1")
        assert 'OTTT' in entry.xref_data
        assert entry.xref_data['OTTT'] == ["OTTHUMT00000397386"]

    def test_it_adds_xref_to_ensembl(self):
        entry = self.entry_from('GENCODE', "ENST00000540226.1")
        assert 'ENST00000540226.1' in entry.xref_data['Ensembl']

    def test_it_gets_all_gencode_entries(self):
        assert len(list(self.data())) == 2898

    def test_it_sets_accession_correctly(self):
        val = self.entries_for("ENST00000540868.1")
        assert len(val) == 2
        assert {(e.accession, e.database) for e in val} == set([
            ('ENST00000540868.1', 'ENSEMBL'),
            ('GENCODE:ENST00000540868.1', 'GENCODE'),
        ])
