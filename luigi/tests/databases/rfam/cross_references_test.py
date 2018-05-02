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

import tempfile
from StringIO import StringIO

import pytest

import databases.rfam.cross_references as cr

from ontologies.data import Term
from ontologies import helpers as ont


@pytest.fixture
def data():
    with open('data/rfam/database_link.tsv', 'r') as raw:
        return list(cr.parse(raw))


def test_can_fetch_and_parse_data(data):
    assert len(data) == 7909


def test_correctly_parses_so_data(data):
    assert data[0] == cr.RfamDatabaseLink(
        'RF00014',
        'SO',
        None,
        'SO:0000378',
        'DsrA_RNA',
    )


def test_correctly_parses_other(data):
    assert data[4] == cr.RfamDatabaseLink(
        'RF00016',
        'snoopy',
        None,
        'Mus_musculus300888;',
        None,
    )


def test_can_fetch_ontology_data(data):
    assert data[0].ontology_term() ==  Term(
        'SO',
        'SO:0000378',
        'DsrA_RNA',
        ('DsrA RNA regulates both transcription, by overcoming '
         'transcriptional silencing by the nucleoid-associated H-NS protein, '
         'and translation, by promoting efficient translation of the stress '
         'sigma factor, RpoS. These two activities of DsrA can be separated '
         'by mutation: the first of three stem-loops of the 85 nucleotide RNA '
         'is necessary for RpoS translation but not for anti-H-NS action, '
         'while the second stem-loop is essential for antisilencing and less '
         'critical for RpoS translation. The third stem-loop, which behaves '
         'as a transcription terminator, can be substituted by the trp '
         'transcription terminator without loss of either DsrA function. The '
         'sequence of the first stem-loop of DsrA is complementary with the '
         'upstream leader portion of RpoS messenger RNA, suggesting that '
         'pairing of DsrA with the RpoS message might be important for '
         'translational regulation.'),
        ['DsrA RNA']
    )


def test_can_extract_all_ontology_terms():
    with open('data/rfam/database_link.tsv', 'r') as raw:
        sample = StringIO()
        for line in raw.readlines()[:10]:
            sample.write(line)
        sample.seek(0)
        assert list(cr.ontology_terms(sample)) == [
            ont.term('SO:0000378'),
            ont.term('SO:0000593'),
            ont.term('GO:0006396'),
            ont.term('GO:0005730'),
        ]
