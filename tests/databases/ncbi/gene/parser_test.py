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

import tempfile

import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import publications as pub

from rnacentral_pipeline.databases.ncbi.gene import fetch
from rnacentral_pipeline.databases.ncbi.gene import helpers
from rnacentral_pipeline.databases.ncbi.gene import parser


def parse(filename):
    with tempfile.NamedTemporaryFile() as tmp:
        with open(filename, 'r') as raw:
            ncrnas = list(helpers.ncrnas(raw))
            fetch.fetch_and_write(ncrnas, tmp)
        tmp.flush()
        tmp.seek(0)
        return list(parser.parse(tmp))


@pytest.fixture(scope='module')
def simple():
    parsed = parse('data/ncbi_gene/simple.txt')
    return {p.primary_id: p for p in parsed}


def test_parses_all_data(simple):
    assert set(simple.keys()) == {
        'NCBI_GENE:113344155',
        'NCBI_GENE:114058594',
        'NCBI_GENE:107303968',
        'NCBI_GENE:112780674',
        'NCBI_GENE:112442518',
        'NCBI_GENE:115038462',
        'NCBI_GENE:107303968',
        'NCBI_GENE:113275511',
        'NCBI_GENE:34399170',
        'NCBI_GENE:108869573',
        'NCBI_GENE:31217049',
        'NCBI_GENE:38465973',
    }


def test_produces_correct_data(simple):
    assert simple['NCBI_GENE:113344155'] == data.Entry(
        primary_id='NCBI_GENE:113344155',
        accession='NCBI_GENE:113344155',
        ncbi_tax_id=3469,
        database='NCBI_GENE',
        sequence='TTTTGGCAGTGATGACTTATACACAGCTTTATGCCAGTTCTGCTTCAGAAAATTCTTGATAGAGAGCTGAAATATATCTGAGTCTAAT',
        regions=[],
        rna_type='snoRNA',
        url='https://www.ncbi.nlm.nih.gov/gene/113344155',
        seq_version='1',
        species='Papaver somniferum',
        common_name='opium poppy',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; '
            'Ranunculales; Papaveraceae; Papaveroideae; Papaver; Papaver '
            'somniferum'
        ),
        xref_data={},
        product="small nucleolar RNA Z223",
        gene='LOC113344155',
        references=[pub.reference(25355515)],
        description='Papaver somniferum (opium poppy) LOC113344155',
    )


@pytest.mark.parametrize('gene_id,rna_type', [
    ('NCBI_GENE:113344155', 'snoRNA'),
    ('NCBI_GENE:107303968', 'lncRNA'),
    ('NCBI_GENE:114058594', 'lncRNA'),
    ('NCBI_GENE:112780674', 'snoRNA'),
    ('NCBI_GENE:112442518', 'tRNA'),
    ('NCBI_GENE:115038462', 'tRNA'),
    ('NCBI_GENE:113275511', 'lncRNA'),
    ('NCBI_GENE:34399170', 'tRNA'),
    ('NCBI_GENE:108869573', 'tRNA'),
    ('NCBI_GENE:31217049', 'tRNA'),
    ('NCBI_GENE:38465973', 'tRNA'),
])
def test_assigns_correct_rna_types(simple, gene_id, rna_type):
    assert simple[gene_id].rna_type == rna_type


# @pytest.mark.parametrize('gene_id,product', [
#     ('NCBI_GENE:112442518', 'tRNA-Leu'),
#     ('NCBI_GENE:115038462', 'tRNA-Glu'),
#     ('NCBI_GENE:34399170', None),
#     ('NCBI_GENE:108869573', 'tRNA-Gly'),
#     ('NCBI_GENE:31217049', 'tRNA-Val'),
#     ('NCBI_GENE:38465973', 'tRNA-Pro')
# ])
# def test_gets_correct_trna_product(simple, gene_id, product):
#     assert simple[gene_id].product == product
