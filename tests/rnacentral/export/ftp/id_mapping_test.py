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

import os
import csv
import tempfile

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import id_mapping as ids
from tests.helpers import run_with_upi_taxid_constraint


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
    assert ids.gene({'gene': gene, 'database': 'PDBE'}) == expected


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


@pytest.mark.parametrize('rna_id,expected', [
    ('URS000018F875_9606', [
        ['URS000018F875', 'ENSEMBL', 'ENST00000523510', 9606, 'lncRNA', 'ENSG00000254166.2'],
        ['URS000018F875', 'GENCODE', 'ENST00000523510', 9606, 'lncRNA', 'CASC19'],
        ['URS000018F875', 'HGNC', 'HGNC:45089', 9606, 'lncRNA', 'PCAT2'],
        ['URS000018F875', 'LNCIPEDIA', 'lnc-FAM84B-15:108', 9606, 'lncRNA', 'lnc-FAM84B-15'],
        ['URS000018F875', 'NONCODE', 'NONHSAT129016.2', 9606, 'lncRNA', 'NONHSAG051247.2'],
        ['URS000018F875', 'REFSEQ', 'NR_119373', 9606, 'lncRNA', 'PCAT2'],
        # ['URS000018F875', 'GENCODE', 'ENST00000523510', 9606, 'lncRNA', 'ENSG00000254166.2'],
    ]),
])
def test_can_create_expected_exports(rna_id, expected):
    entries = run_with_upi_taxid_constraint(
        rna_id,
        'files/ftp-export/id-mapping/query.sql',
        take_all=True
    )
    entries = sorted(ids.as_entry(e) for e in entries)
    assert entries == expected


@pytest.mark.skip()
def test_can_split_by_database():
    def entries():
        rows = [
            ['URS0000000001', 'ENA', 'GU786683.1:1..200:rRNA', '77133', 'rRNA', ''],
            ['URS0000000001', 'ENA', 'GU786868.1:1..200:rRNA', '77133', 'rRNA', ''],
            ['URS0000000001', 'ENA', 'GU786889.1:1..200:rRNA', '77133', 'rRNA', ''],
            ['URS0000000001', 'ENA', 'GU790934.1:1..200:rRNA', '77133', 'rRNA', ''],
            ['URS0000000001', 'ENA', 'GU792481.1:1..200:rRNA', '77133', 'rRNA', ''],
        ]
        for row in rows:
            yield row
            modified = list(row)
            modified[1] = 'new_DB'
            yield modified

    with tempfile.NamedTemporaryFile(suffix='.tsv') as raw:
        base = os.path.dirname(raw.name)
        writer = csv.writer(raw, delimiter='\t')
        writer.writerows(entries())
        raw.flush()

        with open(raw.name, 'r') as r2:
            ids.split_by_database(r2, base)

        ena_file = os.path.join(base, 'ena.tsv')
        db_file = os.path.join(base, 'new_db.tsv')
        assert os.path.exists(ena_file)
        assert os.path.exists(db_file)

        with open(ena_file, 'r') as cur:
            reader = csv.reader(cur, delimiter='\t')
            assert list(reader) == [
                ['URS0000000001', 'ENA', 'GU786683.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'ENA', 'GU786868.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'ENA', 'GU786889.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'ENA', 'GU790934.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'ENA', 'GU792481.1:1..200:rRNA', '77133', 'rRNA', ''],
            ]

        with open(db_file, 'r') as cur:
            reader = csv.reader(cur, delimiter='\t')
            assert list(reader) == [
                ['URS0000000001', 'new_DB', 'GU786683.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'new_DB', 'GU786868.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'new_DB', 'GU786889.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'new_DB', 'GU790934.1:1..200:rRNA', '77133', 'rRNA', ''],
                ['URS0000000001', 'new_DB', 'GU792481.1:1..200:rRNA', '77133', 'rRNA', ''],
            ]

        os.remove(ena_file)
        os.remove(db_file)
