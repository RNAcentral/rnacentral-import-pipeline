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

from rnacentral_pipeline.databases import data as dat
from rnacentral_pipeline.databases.refseq import parser


@pytest.mark.parametrize('filename,count', [
    ('data/refseq/biomol_ncRNA_RNase_MRP_RNA.gbff', 4),
    ('data/refseq/mir-with-several-locations.gbff', 2),
])
def test_can_parse_refseq_files(filename, count):
    with open(filename, 'r') as raw:
        assert len(list(parser.parse(raw))) == count


@pytest.mark.skip()
def test_extracts_correct_chromosome_if_several():
    with open('data/refseq/mir-with-several-locations.gbff', 'r') as raw:
        data = next(parser.parse(raw))
    assert data.chromosome == 'X'


@pytest.mark.xfail()
def test_extracts_correct_coordinates():
    with open('data/refseq/mir-with-several-locations.embl', 'r') as raw:
        data = next(parser.parse(raw))
    assert data.exons == [
        dat.Exon(chromosome='X',
                 primary_start=2609191,
                 primary_end=2609254,
                 complement=False)
    ]


def test_it_can_build_correct_entry():
    with open('data/refseq/biomol_ncRNA_RNase_MRP_RNA.gbff', 'rb') as raw:
        data = list(parser.parse(raw))

    assert len(data) == 4
    data = data[0]
    data = attr.asdict(data)
    assert len(data['references']) == 10
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
            'endoribonuclease (RMRP), RNase MRP RNA'
        ),
        xref_data={
            'GeneID': ['6023'],
            'HGNC': ['HGNC:10031'],
            'MIM': ['157660'],
        },
        note_data={},
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
        keywords='RefSeq',
        optional_id="GeneID:6023",
        product='RNA component of mitochondrial RNA processing endoribonuclease',
        mol_type="transcribed RNA",
    ))


def test_can_build_correct_entries_when_multiple_present():
    with open('data/refseq/mir-with-several-locations.gbff', 'r') as raw:
        data = [attr.asdict(e) for e in parser.parse(raw)]

    assert len(data) == 2
    assert len(data[0]['references']) == 2
    assert len(data[1]['references']) == 2
    assert data[0]['references'] == data[1]['references']

    data[0]['references'] = []
    data[1]['references'] = []

    assert data[0] == attr.asdict(dat.Entry(
        primary_id='NR_106737',
        accession='NR_106737.1:1..64:precursor_RNA',
        ncbi_tax_id=9606,
        database='REFSEQ',
        sequence=(
            'CCCCGGGCCCGGCGTTCCCTCCCCTTCCGTGCGCCAGTGGAGGCCGGGGTGGGGCGGGGCGGGG'
        ),
        exons=[],
        rna_type='precursor_RNA',
        url='https://www.ncbi.nlm.nih.gov/nuccore/NR_106737.1',
        seq_version='1',
        parent_accession='NR_106737',
        is_composite='N',
        description='Homo sapiens microRNA 6089 (MIR6089), precursor RNA',
        xref_data={
            'GeneID': ['102464837'],
            'HGNC': ['HGNC:50179'],
            'miRBase': ['MI0020366'],
        },
        note_data={},
        chromosome='X',
        species="Homo sapiens",
        common_name=u'human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        gene="MIR6089",
        gene_synonyms=[
            'hsa-mir-6089-1',
            'hsa-mir-6089-2',
            'MIR6089-1',
            'MIR6089-2',
        ],
        keywords='RefSeq',
        optional_id="GeneID:102464837",
        product='microRNA 6089',
        mol_type="transcribed RNA",
    ))

    assert data[1] == attr.asdict(dat.Entry(
        primary_id='NR_106737',
        accession='NR_106737.1:39..62:ncRNA',
        ncbi_tax_id=9606,
        database='REFSEQ',
        sequence='GGAGGCCGGGGTGGGGCGGGGCGG',
        exons=[],
        rna_type='miRNA',
        url='https://www.ncbi.nlm.nih.gov/nuccore/NR_106737.1',
        seq_version='1',
        parent_accession='NR_106737',
        is_composite='N',
        description='Homo sapiens microRNA 6089 (MIR6089), miRNA',
        xref_data={
            'GeneID': ['102464837'],
            'HGNC': ['HGNC:50179'],
            'miRBase': ['MI0020366', 'MIMAT0023714'],
        },
        note_data={},
        chromosome='X',
        species="Homo sapiens",
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; '
            'Euteleostomi; Mammalia; Eutheria; Euarchontoglires; '
            'Primates; Haplorrhini; Catarrhini; Hominidae; Homo; '
            'Homo sapiens'
        ),
        gene="MIR6089",
        gene_synonyms=[
            'hsa-mir-6089-1',
            'hsa-mir-6089-2',
            'MIR6089-1',
            'MIR6089-2',
        ],
        keywords='RefSeq',
        optional_id="GeneID:102464837",
        product='hsa-miR-6089',
        mol_type="transcribed RNA",
    ))
