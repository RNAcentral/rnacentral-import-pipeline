# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import os

import pytest
import psycopg2

from rnacentral_pipeline.databases.genecards_suite.core import lookup as lk
from rnacentral_pipeline.databases.genecards_suite.core.data import KnownSequence


@pytest.fixture(scope='module')
def conn():
    return psycopg2.connect(os.environ['PGDATABASE'])


def test_lookups_expected_count(conn):
    urs = [
        'URS00008BD1D3_9606',
        'URS00001EE9F1_9606',
        'URS0000B21DD0_9606',
    ]
    assert len(lk.lookup_chunk(conn, urs)) == len(urs)


def test_produces_correct_data(conn):
    urs = ['URS0000D58B85_9606']
    results = lk.lookup_chunk(conn, urs)
    assert len(results) == 1
    assert results[0] == KnownSequence(
        rna_id='URS0000D58B85_9606',
        rna_type='lncRNA',
        sequence=(
            'GCCCCCGCCCATCAGCCCAGCCTGCAAGGAGGCGCCACCGGACCTGCA'
            'CACTCAGGAGCAGACTCAGCTAGATTGGAAGAGAGGATGGGAAGACAATGGATGACCTAC'
            'CTGTGGAGAGAAACTACCCATTCCAGAACCTCCTCTCTGCTGAAAGCTGAACACTCAATG'
            'GGAAGACCTGTTTGCAGAGAGGAGCCAGCCACTCCAGGGTCTCCGTTCTGCTAGAAGCTG'
            'AACACTCGTTGGGACACCCTTGCTGTGGAAAGGACATCCTTACTGTGGAAATGAGCTACT'
            'AACTGCAGGAGCTGAGGAGCTGCACAGCAAGGTTAAAAGAGCTGTAACACAAACAGGGCT'
            'GCAACATGCCCCTTGCTCCCCACAGGGAGAGAAGAGCTCTGGCCCTCGGAGAAGCCCAGA'
            'CCTGGGAGCTCCTTGAGCCCGGGCTGTGACTCCCTCTTTGGGGCCCTGGTTGGCGTCACT'
            'GCATTCGCCAGTGCCACTGTTGGAAGCTGCTTGTGATGCGCCTGGTCCAGGGGGAAGCTG'
            'TTTGTTGTGTGCCTGGTCCAGCCACCTCATGGAGAGCCTGTGCTGGCACCTGGAGCTGCC'
            'CAACCTGGGCAGCAGCTTGTGTGTGTGACTGCACAGTGGCCACGCTTGCTCACACACCCC'
            'TAGCACACCCCTCTGCTCCACCCGTCTCAATCTCCCTTGGATCCGGGATCCAGTCCTTAC'
            'ATGAAGCCTACTCAGACCGACATCACTGTATCATTTTTTTGTATTCTGAAGTCCTTCAGC'
            'ACTTGATTCCTTTAAGTTTAATCTTGGTATTCTATTGGGAATCATATGGAAATGTACAAT'
            'GATATTTATTGAATGAACAGGAAACAGGGAAAGACATTGACCAGAAAAAGTGTTTATTC'
        ),
        description='Homo sapiens (human) non-protein coding lnc-KLRG1-9:4',
    )


def test_can_lookup_and_index(conn):
    urs = [
        'URS00008BD1D3_9606',
        'URS00001EE9F1_9606',
        'URS0000B21DD0_9606',
    ]
    data = lk.lookup(urs, conn)
    assert len(data) == len(urs)
    for rna_id in urs:
        assert rna_id in data
