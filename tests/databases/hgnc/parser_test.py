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

import os
import json
from cStringIO import StringIO
from collections import Counter

import attr

import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases import mapping

from rnacentral_pipeline.databases.hgnc import parser as hgnc
from rnacentral_pipeline.databases.hgnc import helpers


def on_selected(index, handler):
    with open('data/hgnc/data.json') as raw:
        data = json.load(raw)
        docs = data['response']['docs']
        if isinstance(index, int):
            item = data['response']['docs'][index]
        elif isinstance(index, basestring):
            item = next(d for d in docs if d['hgnc_id'] == index)

        assert item
        data['response']['docs'] = [item]

    selected = StringIO(json.dumps(data))
    return handler(selected)


@pytest.fixture
def partials():
    with open('data/hgnc/data.json') as raw:
        return hgnc.parse_partials(raw)


def partial(index):
    return next(on_selected(index, hgnc.parse_partials))


@pytest.fixture
def entries():
    with open('data/hgnc/data.json') as raw:
        return hgnc.parse(raw, os.environ.get('PGDATABASE'))


def entry(index):
    def handler(raw):
        return hgnc.parse(raw, os.environ.get('PGDATABASE'))
    return list(on_selected(index, handler))


@pytest.mark.slowtest
def test_can_create_a_partial_parse_of_all_data(partials):
    assert len(list(partials)) == 7257


def test_can_correctly_partially_parse_an_entry():
    assert attr.asdict(partial(0)) == attr.asdict(mapping.SequencelessEntry(
        primary_id="A1BG-AS1",
        accession="HGNC:37133",
        ncbi_tax_id=9606,
        database='HGNC',
        exons=[],
        rna_type='lncRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:37133',
        seq_version='1',
        xref_data={
            'RefSeq': ["NR_015380"],
            'UCSC': ["uc002qse.3"],
            'LNCipedia': ["A1BG-AS1"],
            'Ensembl': ["ENSG00000268895"],
            'ENA': ["BC040926"],
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene="A1BG-AS1",
        description="Homo sapiens (human) A1BG antisense RNA 1 (A1BG-AS1)",
        gene_synonyms=[
            "NCRNA00181",
            "A1BGAS",
            "A1BG-AS",
        ],
        references=[],
    ))


@pytest.mark.slowtest
def test_it_can_correctly_parse_and_map_all_entries(entries):
    assert len(list(entries)) == 7248


def test_can_correctly_map_an_entry():
    parsed = entry(0)
    assert len(parsed) == 1
    assert attr.asdict(parsed[0]) == attr.asdict(data.Entry(
        primary_id="A1BG-AS1",
        accession="HGNC:37133",
        ncbi_tax_id=9606,
        database='HGNC',
        sequence=(
            'ATTTTTAGTAGAGACGGGGTTTCGTCATGTTGGGTAGACTGGTCTCGATCTCTTGACCTC'
            'ATGATCCGCTGGCCTCAGCCTCCCAAAGTGCTGGGATTATAGGCGTGAGCCACCGCACCC'
            'GGCCTCACCCTTCCATTCTTTGGGCATCCTCTGCCCAAACTCCTTTGCAACCATCAGACC'
            'CAGGCAACCGTGTGGCCAGCTTCCTGGACCTGGTCCTGGGCAAGGACATTCAGCATGGAG'
            'GTCAGGATGAATGTGGTCATGCCCCCTGTATCCTCCTCCCCCAGCAGCCCCCAGAGATGG'
            'TTCCCACAGTCACCGAGCTCCTCAATGGTCACAGTAGCGCTGGGCTCAGAGAGGGCGCCT'
            'TCCCCATCGGTCCGGTAGCTGCAGCTGTAGTTGCCAGGCTGATGGACTGGAAAGGTGGCC'
            'TCCACATCCTCCTGGGCCTCAGGCACCTCCAGAAACTCATGGTCGCCCTCCCGCCTCAGC'
            'AGAAAAGTCACACCCCGCAGCACACCTCGGCACACTGCTGTTGTTTTCAGGCCGGGGGTG'
            'ATCCAGGACACTGGCGCCATCGAGAGCCAGGGAGCAGGCAAGGACTCTGTGGGTGGGGTG'
            'GAGTTGAAGAGTGCTTCTGCATAACAGGAGACGTGGGAAGCCCAGCAGAACCTGGGCCCC'
            'ACACTCATAAGACTGTGAAGGAGACAGGGATGGAGGGCATTGCTGGGAAAAGGCAGGCAG'
            'GCACCAGTGGAAACAAGGCATACGGCAAGAGAAAAAGAGTGTGTGGGAGCTGTACGGCAT'
            'GCCGGGCCTGCTGGGGACACAACCATGTATTCCTTACAACACAGACATCATCCCCATTTG'
            'ACAGATGGGAAGATTGGGGCTTGGGGCTCAGAGACATCGGGTGACTTGGAGGAAGGGAGA'
            'CCCCCCAGTCAACAACCTACCCACTGCCTCCTGGCCCCCAATTCATGCCCCGGCTCACTT'
            'GGCCCTGTCAGCTCCAGGAGCTTGCTCAGCTGGGTCCATCCTGTGGACAAGCCCGAGCGG'
            'CAGCGGTAGCGGCCCTGGGTGTCACCCGTCAGCAGGAACTGGTGCTTGATGGCAGGGACC'
            'CAGACCCCAGGAGGCCTCACCCCACAGCAAGAGAAAGACCACGAGCATGGACATGATGGT'
            'CGCGCTCACTCCGGTGCAGTGAGTGTCTGGGTGATGGATCCTAGAACCTGTTCTAGAAGG'
            'AGACGCTGAGTTCTCCCCACTCTACCTGGAGCTGCTGAGGTCAAGTGGAGAGTGCCGGGG'
            'CTGAAGCCCAGCGAGCCCCTCCCTGTCACACCACACTGCACATGTGACACCAAAGGCATT'
            'TATTGAGCATCTACTGTGTGCCAGATGCAAAGAAGACAGTGGGATAACTGGGCAAGGAGC'
            'ATTCCAGGCAGAGGGAACAGCCAGTGCAAACATGGGAGATGCGGGTGCAGTGAGAAGGAG'
            'CAGAGGAGACCAGGAGTTGAGGCGACTGCAGAGCGATTGCAGAGCTCATAGACGTGATGA'
            'AGACCAGCTGGGCACTGCAGCGGCTCTGGCTTCTGCTGTGAGGGACGGGGAGCTGTGGAG'
            'ACTCCCTCTGTTTGGTCCCATTTCCATATTGCTGAGCACCTGCTGCATGCTCTCTGTCCT'
            'GCTCAGAGCCAGTACCTGGATGGAGGCTGTGTGTTCTGGGTGGACTGGGGGAGAATGCAT'
            'TCTTTGCAGAAGAGACAACCAGGTGCTGAGGCCTGAAGTCATAACCAGACCTGGGGCACT'
            'CAGACAAGGACTGGTGGCCAGGCCAGAGGAGCAGAGCTCAGGGTGTGCACAGAACTCAGA'
            'AGTCAGACCACTGAACAGTGACAGAACCTTTCAGCCCATTGGCAATGAGGCTGCCCAGGC'
            'AGCCAGAGGCCTGGTCAAATCTGAGGTCTGCAGAGATGGAGCTGTCATTCTGTGTTTCCT'
            'TTGGCAGTCCCAACACCAGCCTCGTTGCACCCTCCTCTTAGCTTCCCTGGGTTCTCCAGC'
            'CCTAAGGGTGGTAGCAGCTTCCTGCAAGTATCCAGCTCTGAGATTCTGCAACATCCACTT'
            'TTGCTCTCTGAGCCTTGCAAAACCTGCACAATCAGTCCCTAACTTATATCCCCTCTGTTT'
            'GAAATACCTAGTGTGGTTTCTATTTCCTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
            'AAAAAAAAAAAAAAA'
        ),
        exons=[],
        rna_type='lncRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:37133',
        seq_version='1',
        xref_data={
            'RefSeq': ["NR_015380"],
            'UCSC': ["uc002qse.3"],
            'LNCipedia': ["A1BG-AS1"],
            'Ensembl': ["ENSG00000268895"],
            'ENA': ["BC040926"],
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene="A1BG-AS1",
        description="Homo sapiens (human) A1BG antisense RNA 1 (A1BG-AS1)",
        gene_synonyms=[
            "NCRNA00181",
            "A1BGAS",
            "A1BG-AS",
        ],
        references=[],
    ))


def test_knows_when_sequence_does_not_map():
    assert not entry(9)  # ABCF1-DT


def test_can_correctly_partially_parse_a_tRNA(partials):
    parsed = entry('HGNC:34987')
    assert len(parsed) == 1
    assert attr.asdict(parsed[0]) == attr.asdict(data.Entry(
        primary_id='TRE-TTC1-1',
        accession='HGNC:34987',
        ncbi_tax_id=9606,
        database='HGNC',
        sequence=(
            'TCCCATATGGTCTAGCGGTTAGGATTCCTGGTTTTCACCCAGGTGGCCCGGGTTCGACTC'
            'CCGGTATGGGAA'
        ),
        exons=[],
        rna_type='tRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:34987',
        seq_version='1',
        xref_data={
            'ENA': ['HG983778'],
            'GtRNAdb': ['tRNA-Glu-TTC-1-1'],
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene='TRE-TTC1-1',
        description='Homo sapiens (human) transfer RNA-Glu (TTC) 1-1 (TRE-TTC1-1)',
        gene_synonyms=['TRNAE26', 'TRE-TTC1-2'],
        references=[],
    ))


@pytest.mark.skip
def test_can_correctly_partially_parse_a_mito_tRNA(partials):
    assert False


def test_can_map_to_ensembl():
    parsed = entry('HGNC:41309')
    assert len(parsed) == 1
    assert attr.asdict(parsed[0]) == attr.asdict(data.Entry(
        primary_id='ATP2B2-IT1',
        accession='HGNC:41309',
        ncbi_tax_id=9606,
        database='HGNC',
        sequence=(
            'CTGCAAATCCCTCCCACTCACCAAACTGCTGCCGCCTGCAATTCCACAATGAGCACTTTA'
            'GTCATCAAATCATAAAATGTCAGCTTGTAATTAGCCTTGGAACCTGTCTAGTTCAATTTC'
            'CCAGCTCAGAGTGCTCTTGGAAGATACCTGAACGGGCCCTTTATTTTAGAGTGAATTAGT'
            'GGCAGAGCTGGAACTCGAGTCTGACCCTGCACAGTGAAGGACGCTAGAGGCATCTGCAGC'
            'TTCTTCAATGGGAAGTTGGGGCATCATATGCCTTTTTGAAATTATGTGTATTGGTAGGTT'
            'ATATGATGTGTGAATTTTCTTGGAAGATAAAGAAAAAAAGTAAAATGATGTTGGATCCAA'
            'TGGGACCAAGAGCCACCAGCACAGTGGATAGAAACAAAGGCCCTGTACTGAATCCTGCCC'
            'TGCCTCCTAATTAGCTAGGACCTTGCATATGAGCATGAACATCAGCTCTTAC'
        ),
        exons=[],
        rna_type='lncRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:41309',
        seq_version='1',
        xref_data={
            'Ensembl': ['ENSG00000236999'],
            'LNCipedia': ['ATP2B2-IT1'],
            'UCSC': ['uc062gsp.1'],
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene='ATP2B2-IT1',
        description='Homo sapiens (human) ATP2B2 intronic transcript 1 (ATP2B2-IT1)',
        gene_synonyms=[],
        references=[]
    ))

def test_can_correctly_partially_parse_a_sequence_in_mitochondira(partials):
    parsed = entry('HGNC:35035')
    assert len(parsed) == 1
    assert attr.asdict(parsed[0]) == attr.asdict(data.Entry(
        primary_id='NMTRP-TGG1-1',
        accession='HGNC:35035',
        ncbi_tax_id=9606,
        database='HGNC',
        sequence=(
            'TAGGACTTGGTGTAATAGGTAGCACGAAGAGATTTGGATTCTCAGGGGTAGGTTCAATTC'
            'CTATAGTTCTGGAA'
        ),
        exons=[],
        rna_type='tRNA',
        url='https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=HGNC:35035',
        seq_version='1',
        xref_data={
            'ENA': [u'HG984168'],
            'RefSeq': [u'NG_008520'],
            'GtRNAdb': ['nmt-tRNA-Pro-TGG-1-1']
        },
        species='Homo sapiens',
        common_name='human',
        lineage=(
            "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; "
            "Euteleostomi; Mammalia; Eutheria; Euarchontoglires; "
            "Primates; Haplorrhini; Catarrhini; Hominidae; Homo; "
            "Homo sapiens"
        ),
        gene='NMTRP-TGG1-1',
        description='Homo sapiens (human) nuclear-encoded mitochondrial transfer RNA-Pro (TGG) 1-1 (NMTRP-TGG1-1)',
        gene_synonyms=['TRNAP23P'],
        references=[],
    ))


@pytest.mark.slowtest
def test_can_correctly_get_all_mito_sequences(partials):
    mito = [p for p in partials if p.organelle == 'Mitochondrion']
    assert len(mito) == 25


@pytest.mark.skip
def test_can_correctly_load_with_references(partials):
    assert False


@pytest.mark.slowtest
def test_fetches_correct_rna_counts(partials):
    counts = Counter(e.rna_type for e in partials)
    assert dict(counts) == {
        'Y_RNA': 4,
        'other': 124,
        'lncRNA': 3913,
        'precursor_RNA': 1878,
        'misc_RNA': 30,
        'rRNA': 60,
        'scRNA': 3,
        'snRNA': 37,
        'snoRNA': 567,
        'tRNA': 637,
        'vault_RNA': 4,
    }
