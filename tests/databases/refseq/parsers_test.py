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

from databases import data as dat
from databases.refseq import parsers


@pytest.mark.parametrize('filename,count', [
    ('data/refseq/biomol_ncRNA_RNase_MRP_RNA.dat', 4),
    ('data/refseq/mir-with-several-locations.embl', 1),
])
def test_can_parse_refseq_files(filename, count):
    with open(filename, 'r') as raw:
        assert len(list(parsers.parse(raw))) == count


def test_extracts_correct_chromosome_if_several():
    with open('data/refseq/mir-with-several-locations.embl', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.chromosome == 'X'


@pytest.mark.xfail()
def test_extracts_correct_coordinates():
    with open('data/refseq/mir-with-several-locations.embl', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.exons == [
        dat.Exon(chromosome='X',
                 primary_start=2609191,
                 primary_end=2609254,
                 complement=False)
    ]


def test_can_correctly_assign_refseq_db():
    with open('data/refseq/biomol_ncRNA_RNase_MRP_RNA.dat', 'r') as raw:
        data = next(parsers.parse(raw))
    assert data.database == 'REFSEQ'


def test_can_correctly_assign_external_id():
    with open('data/test_refseq_product.ncr', 'rb') as raw:
        data = next(parsers.parse(raw))
    assert data.primary_id == 'NR_029991'


def test_can_correctly_assign_optional_id():
    with open('data/test_refseq_product.ncr', 'rb') as raw:
        data = next(parsers.parse(raw))
    assert data.optional_id == "GeneID:100033549"


def test_it_can_build_correct_entry():
    with open('data/refseq/biomol_ncRNA_RNase_MRP_RNA.dat', 'rb') as raw:
        data = next(parsers.parse(raw))

    data = attr.asdict(data)
    assert len(data['references']) == 12
    data['references'] = []

    assert data == attr.asdict(dat.Entry(
        primary_id='NR_003051',
        accession='NR_003051.3:1..277:ncRNA',
        ncbi_tax_id=9606,
        database='REFSEQ',
        sequence=(
            'GGTTCGTGCTGAAGGCCTGTATCCTAGGCTACACACTGAGGACTCTGTTCCTCCCCTTTC'
            'CGCCTAGGGGAAAGTCCCCGGACCTCGGGCAGAGAGTGCCACGTGCATACGCACGTAGAC'
            'ATTCCCCGCTTCCCACTCCAAAGTCCGCCAAGAAGCGTATCCCGCTGAGCGGCGTGGCGC'
            'GGGGGCGTCATCCGTCAGCTCCCTCTAGTTACGCAGGCAGTGCGTGTCCGCGCACCAACC'
            'ACACGGGGCTCATTCTCAGCGCGGCTGTAAAAAAAAA'
        ),
        exons=[],
        rna_type='RNase_MRP_RNA',
        url='https://www.ncbi.nlm.nih.gov/nuccore/NR_003051.3',
        seq_version='3',
        parent_accession='NR_003051',
        is_composite='N',
        description=(
            'Homo sapiens RNA component of mitochondrial RNA processing '
            'endoribonuclease (RMRP), RNase MRP RNA.'
        ),
        project='PRJEB6684',
        xref_data={
            'GeneID': ['6023'],
            'HGNC': ['HGNC:10031'],
            'MIM': ['157660'],
            'ena_refs': {
                'REFSEQ': ('NR_003051', None)
            }
        },
        note_data={
            'ontology': [
                'ECO:0000345',
                'GO:0000171',
                'GO:0000172',
                'SO:0000385',
            ],
        },
        chromosome='9',
        species="Homo sapiens",
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        gene="RMRP",
        gene_synonyms=["CHH", "NME1", "RMRPR", "RRP2"],
        keywords='RefSeq; RNAcentral; TPA; TPA:specialist_db',
        optional_id="GeneID:6023",
        product='RNA component of mitochondrial RNA processing endoribonuclease',
        mol_type="transcribed RNA",
    ))
