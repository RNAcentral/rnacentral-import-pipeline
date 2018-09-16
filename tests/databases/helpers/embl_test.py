# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.helpers import embl


def test_gets_the_taxid(human_12):
    assert embl.taxid(self.record) == 9606


def test_it_can_get_gene(human_12):
    assert embl.gene(self.features['misc_RNA', 'ENSG00000256263.1']) == \
        'ENSG00000256263.1'


def test_it_can_get_the_locus_tag(human_12):
    assert embl.locus_tag(self.features['gene', 'ENSG00000256948.1']) == \
        'RP11-598F7.3'


def test_it_knows_if_something_is_a_gene(human_12):
    assert embl.is_gene(self.features['gene', 'ENSG00000256948.1']) is True
    assert embl.is_gene(self.features['misc_RNA', 'ENSG00000120645.11']) is \
        False


def test_it_can_get_grouped_xref_data(human_12):
    feature = self.features['misc_RNA', 'ENSG00000256948.1']
    assert embl.xref_data(feature) == {
        "Vega_transcript": ["OTTHUMT00000397386"],
        "UCSC": ["uc058jnh.1"],
        "Clone_based_vega_transcript": ["RP11-598F7.3-001"],
        "OTTT": ["OTTHUMT00000397386"],
        "RNAcentral": ["URS000010A1F5"],
    }
