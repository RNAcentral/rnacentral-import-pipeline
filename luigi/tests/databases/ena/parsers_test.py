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

import attr

import pytest

from databases import data
from databases.ena.parsers import parse


@pytest.fixture
def parsed_data():
    with open('data/wgs_aacd01_fun.ncr', 'rb') as raw:
        return list(parse(raw))


def test_creates_all_entries(parsed_data):
    assert len(list(parsed_data)) == 188


def test_creates_simple_entry(parsed_data):
    assert attr.asdict(parsed_data[0]) == attr.asdict(data.Entry(
        primary_id='',
        accession='AACD01000002.1:101667..101773:tRNA',
        ncbi_tax_id=227321,
        database='ENA',
        sequence='GCCCGGATGGTCTAGTGGTATGATTCTCCCTTCGGGGGCAGCGCCCGGTACATAATAACATGTATCAGAAATGGGAGAGGTCCCGCGTTCGAATCGCGGTTCGGGCC',
        exons=[
        ],
        rna_type='tRNA',
        url='',
        seq_version='1',
        chromosome='VIII',
        species='Aspergillus nidulans FGSC A4',
        lineage=(
            'Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; '
            'Eurotiomycetes; Eurotiomycetidae; Eurotiales; Aspergillaceae; '
            'Aspergillus; Aspergillus nidulans FGSC A4'
        ),
        product='tRNA-Pro',
        parent_accession='AACD01000002',
        project='PRJNA130',
        keywords='WGS',
        division='FUN',
        description='Aspergillus nidulans FGSC A4 tRNA-Pro',
        mol_type='genomic DNA',
        is_composite='N',
        references=[

        ],
    ))
