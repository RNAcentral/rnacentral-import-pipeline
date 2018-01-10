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

from internal.export.ftp import id_mapping as ids

from tasks.config import db
from tests.tasks.export.ftp.fasta.active_test import count


def test_can_create_accession_for_pdb():
    data = {
        'database': 'PDBE',
        'external_id': '1S72',
        'optional_id': '1',
        'accession': '1S72_1',
    }
    assert ids.accession(data) == '1S72_1'


def test_can_create_accession_for_hgnc():
    data = {
        'database': 'HGNC',
        'external_id': 'A',
        'optional_id': '1',
        'accession': 'HGNC:A',
    }
    assert ids.accession(data) == 'HGNC:A'


def test_can_create_accession_for_ena():
    data = {
        'database': 'ENA',
        'external_id': None,
        'optional_id': None,
        'accession': 'EU334494.1:47..232:misc_RNA',
    }
    assert ids.accession(data) == 'EU334494.1:47..232:misc_RNA'


def test_can_create_generic_accession():
    data = {
        'database': 'ENSEMBL',
        'external_id': 'ENSSSCT00000052986',
        'optional_id': 'ENSSSCG00000037543.1',
        'accession': 'ENSSSCT00000052986.1',
    }
    assert ids.accession(data) == 'ENSSSCT00000052986'


@pytest.mark.parametrize('gene,expected', [
    ('A', 'A'),
    ('something or other', 'something or other'),
    ('first\tsecond', 'first second'),
])
def test_can_generate_gene(gene, expected):
    assert ids.gene({'gene': gene}) == expected


@pytest.mark.parametrize('name,expected', [
    ('PDBE', 'PDB'),
    ('HGNC', 'HGNC')
])
def test_can_create_databases(name, expected):
    assert ids.database({'database': name}) == expected


def test_as_entry_works_correctly():
    raw = {
        'upi': 'a',
        'database': 'PDBE',
        'external_id': '1S72',
        'optional_id': '1',
        'taxid': 102,
        'rna_type': 'rRNA',
        'gene': None,
    }
    assert ids.as_entry(raw) == [
        'a',
        'PDB',
        '1S72_1',
        102,
        'rRNA',
        '',
    ]


@pytest.mark.slowtest
def test_complete_sql_produces_correct_numbers():
    assert count(ids.COMPLETE_SQL) == 28306193


@pytest.mark.slowtest
def test_complete_produces_correct_count():
    assert len(list(ids.complete(db()))) == 28306193


@pytest.mark.slowtest
def test_example_produces_correct_data():
    from pprint import pprint
    example = list(ids.example(db()))
    pprint(example)
    assert example == [
        ['URS0000000001', 'ENA', 'GU786683.1:1..200:rRNA', '77133', 'rRNA', ''],
        ['URS0000000001', 'ENA', 'GU786868.1:1..200:rRNA', '77133', 'rRNA', ''],
        ['URS0000000001', 'ENA', 'GU786889.1:1..200:rRNA', '77133', 'rRNA', ''],
        ['URS0000000001', 'ENA', 'GU790934.1:1..200:rRNA', '77133', 'rRNA', ''],
        ['URS0000000001', 'ENA', 'GU792481.1:1..200:rRNA', '77133', 'rRNA', '']
    ]
