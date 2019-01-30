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
    ('URS0000672F0E_7955', [
        ['URS0000672F0E', 'ENSEMBL', 'ENSDART00000171022', 7955, 'Y_RNA', 'ENSDARG00000100903.1'],
        # ['URS0000672F0E', 'ENSEMBL', 'ENSDART00000171022.1', 7955, 'Y_RNA', 'ENSDARG00000100903.1'],
        ['URS0000672F0E', 'RFAM', 'RF00019', 7955, 'Y_RNA', ''],
        ['URS0000672F0E', 'RFAM', 'RF00019', 7955, 'Y_RNA', ''],
    ])
])
def test_can_create_expected_exports(rna_id, expected):
    entries = run_with_upi_taxid_constraint(
        rna_id,
        'files/ftp-export/id-mapping/id_mapping.sql',
        take_all=True
    )
    entries = sorted(ids.as_entry(e) for e in entries)
    assert entries == expected
