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

from databases.data import Exon
from databases.data import Entry
from databases.rgd import helpers as rgd


@pytest.fixture
def simple_entry():
    with open('data/rgd/rat_ncrna.tsv', 'rb') as raw:
        return next(rgd.as_rows(raw))


@pytest.fixture
def tricky_entry():
    with open('data/rgd/rat_ncrna.tsv', 'rb') as raw:
        for entry in rgd.as_rows(raw):
            if entry['GENE_RGD_ID'] == '10401674':
                return entry
    raise ValueError("What!")


@pytest.fixture
def rat_ncrna():
    with open('data/rgd/rat_ncrna.tsv', 'rb') as raw:
        return list(rgd.as_rows(raw))


@pytest.fixture
def rat_protein():
    with open('data/rgd/rat_protein.tsv', 'rb') as raw:
        return list(rgd.as_rows(raw))


@pytest.fixture
def sequences():
    with rgd.indexed('data/rgd/sequences.fa.gz') as indexed:
        yield indexed


def test_can_fetch_all_organisms():
    assert rgd.known_organisms() == ['rat']


def test_can_fetch_accession_from_entry(simple_entry):
    assert rgd.accession(simple_entry) == 'RRID:RGD_5687330'


def test_can_fetch_accession_from_entry_with_index(simple_entry):
    assert rgd.accession(simple_entry, 0) == 'RRID:RGD_5687330:0'


def test_can_determine_taxid(simple_entry):
    assert rgd.taxid(simple_entry) == 10116


@pytest.mark.parametrize('entry', rat_ncrna())
def test_can_detect_if_is_ncrna(entry):
    assert rgd.is_ncrna(entry) is True


@pytest.mark.parametrize('entry', rat_protein())
def test_can_detect_if_not_ncrna(entry):
    assert rgd.is_ncrna(entry) is False


def test_can_generate_url(simple_entry):
    assert rgd.url(simple_entry) == 'https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=5687330'


def test_can_generate_exons(simple_entry):
    assert rgd.exons(simple_entry) == [Exon(
        chromosome='chr14',
        primary_start=46643222,
        primary_end=46645093,
        complement=False,
    )]


def test_can_generate_reasonable_description(simple_entry):
    assert rgd.description(simple_entry) == 'Rattus norvegicus 18S ribosomal RNA (Rn18s)'


def test_can_build_xrefs(simple_entry):
    assert rgd.xref_data(simple_entry) == {
        'genbank': ['AABR07015078', 'AH001747', 'DQ623540', 'NR_046237', 'V01270'],
        'ncbi_gene': ['100861533'],
    }


def test_can_fetch_sequence(simple_entry):
    with rgd.indexed('data/rgd/sequences.fa.gz') as indexed:
        assert rgd.sequences_for(simple_entry, indexed) == ['TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGGACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCTCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGTCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACTGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA']


def test_fails_without_existing_sequence(simple_entry, sequences):
    entry = dict(simple_entry)
    entry['GENE_RGD_ID'] = 'something-made-up'
    with pytest.raises(Exception):
        rgd.sequence(entry, sequences)


def test_can_build_correct_entry(simple_entry, sequences):
    assert attr.asdict(rgd.as_entry(simple_entry, sequences)) == attr.asdict(Entry(
        primary_id='5687330',
        accession='RRID:RGD_5687330',
        ncbi_tax_id=10116,
        database='RGD',
        sequence='TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGGACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGCCCTCCGTGGCGGCGACGACCCATTCGAACGTCTGCCCTATCAACTTTCGATGGTAGTCGCCGTGCCTACCATGGTGACCACGGGTGACGGGGAATCAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCAAGGAAGGCAGCAGGCGCGCAAATTACCCACTCCCGACCCGGGGAGGTAGTGACGAAAAATAACAATACAGGACTCTTTCGAGGCCCTGTAATTGGAATGAGTCCACTTTAAATCCTTTAACGAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTTGGGAGCGGGCGGGCGGTCCGCCGCGAGGCGAGTCACCGCCCGTCCCCGCCCCTTGCCTCTCGGCGCCCCCTCGATGCTCTTAGCTGAGTGTCCCGCGGGGCCCGAAGCGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGCCCGAGCCGCCTGGATACCGCAGCTAGGAATAATGGAATAGGACCGCGGTTCTATTTTGTTGGTTTTCGGAACTGAGGCCATGATTAAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACCAGAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCATAAACGATGCCGACTGGCGATGCGGCGGCGTTATTCCCATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACCCGGCCCGGACACGGACAGGATTGACAGATTGATAGCTCTTTCTCGATTCCGTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGATAACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCGAGCGGTCGGCGTCCCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACCCGAGATTGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCTACACTGACTGGCTCAGCGTGTGCCTACCCTACGCCGGCAGGCGCGGGTAACCCGTTGAACCCCATTCGTGATGGGGATCGGGGATTGCAATTATTCCCCATGAACGAGGAATTCCCAGTAAGTGCGGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGCCCTCGGATCGGCCCCGCCGGGGTCGGCCCACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTGACTATCTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA',
        exons=[Exon(
            chromosome='chr14',
            primary_start=46643222,
            primary_end=46645093,
            complement=False,
        )],
        rna_type='rRNA',
        url='https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=5687330',
        seq_version='1',
        xref_data={
            'genbank': ['AABR07015078', 'AH001747', 'DQ623540', 'NR_046237', 'V01270'],
            'ncbi_gene': ['100861533'],
        },
        gene='Rn18s',
        locus_tag='Rn18s',
        description='Rattus norvegicus 18S ribosomal RNA (Rn18s)',
        references=[
        ]
    ))


def test_produces_none_for_no_sequence(tricky_entry, sequences):
    assert rgd.as_entry(tricky_entry, sequences) is None
    # assert attr.asdict(rgd.as_entry(tricky_entry, sequences)) == attr.asdict(Entry(
    #     primary_id='10401674',
    #     accession='RRID:RGD_10401674',
    #     ncbi_tax_id=10116,
    #     database='RGD',
    #     sequence='',
    #     exons=[],
    #     rna_type='ncrna',
    #     url='https://rgd.mcw.edu/rgdweb/report/gene/main.html?id=10401674',
    #     seq_version='1',
    #     xref_data={

    #     },
    #     gene='Carmn',
    #     locus_tag='Carmn',
    #     gene_synonyms=['Mir143hg'],
    #     description='Rattus norvegicus cardiac mesoderm enhancer-associated non-coding RNA (Carmn)',
    # ))


def test_can_create_entries_for_all_ncrna(rat_ncrna, sequences):
    entries = [rgd.as_entry(e, sequences) for e in rat_ncrna]
    assert len(entries) == 15
    assert len([e for e in entries if e]) == 14
