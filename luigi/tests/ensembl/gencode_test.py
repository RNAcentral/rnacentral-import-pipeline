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

import pytest

from ensembl.gencode import Gencode

from tests.ensembl.helpers import Base, FeatureParsingTest, BothParsingTest


class BasicGencodeTest(FeatureParsingTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = Gencode

    def test_knows_if_entry_is_gencode(self):
        feature = self.features[('misc_RNA', 'ENSG00000256263.1')]
        assert self.importer.is_gencode(feature) is True

    def test_knows_if_entry_is_not_gencode(self):
        feature = self.features[('mRNA', 'ENSG00000120645.11')]
        assert self.importer.is_gencode(feature) is False


class OverridingEnsemblDataTest(FeatureParsingTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = Gencode

    def test_it_uses_correct_primary_id(self):
        assert self.importer.gencode_primary_id('CDS', "ENSG00000120645.11") == 'OTTHUMT00000397383'

    @pytest.mark.skip()
    def test_it_creates_valid_accession(self):
        pass

    def test_it_removes_xref_to_gencode(self):
        feature = self.features[('CDS', "ENSG00000120645.11")]
        assert "OTTT:OTTHUMT00000397383" not in self.gencode_xrefs(feature)
        assert "OTTP:OTTHUMP00000237527" not in self.gencode_xrefs(feature)

    def test_it_adds_xref_to_ensembl(self):
        assert 'Ensembl:ENSG00000120645.11' in self.gencode_xrefs('CDS', "ENSG00000120645.11")

    @pytest.mark.skip()
    def test_sets_correct_references(self):
        pass


class CompleteGencodeTest(BothParsingTest):
    filename = 'data/Homo_sapiens.GRCh38.87.chromosome.12.dat'
    importer_class = Gencode

    def test_it_gets_all_gencode_entries(self):
        assert len(list(self.importer.data(self.filename))) == 1000

    def test_it_sets_accession_to_transcript_id(self):
        entries = list(self.rnacentral_entries('misc_RNA', 'ENSG00000255746.1'))
        assert len(entries) == 1
        assert {(e['accession'], e['database']) for e in entries} == set([
            ('ENST00000540868.1', 'ENSEMBL'),
            ('OTTHUMT00000421260', 'GENCODE'),
        ])
