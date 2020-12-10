# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases import data
import rnacentral_pipeline.databases.helpers.publications as pubs

from rnacentral_pipeline.databases.silva import parser


@pytest.mark.parametrize('filename,count', [
    ('data/silva/sample.tsv', 9),
])
def test_parses_all_data(filename, count):
    with open(filename, 'r') as raw:
        assert len(list(parser.parse(raw))) == count


def test_parses_data_correctly():
    with open('data/silva/sample.tsv', 'r') as raw:
        val = next(parser.parse(raw))

    assert val == data.Entry(
        primary_id='SILVA:FN662328.1:1..957',
        accession='SILVA:FN662328.1:1..957',
        ncbi_tax_id=100272,
        database='SILVA',
        sequence=(
            'GGGCGCGTGGGTGACGAAGGTCTTCGGATTGTAAAGCCCTTTTCTGGGGGAAGATGATGA'
            'CGGTACCCCAGGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGG'
            'TGCGAGCGTTGCCCGGATTTATTGGGCGTAAAGGGCGCGTAGGCGGTCACGCACGTCCGT'
            'TGTGAAATCGCTCGGCTCAACTGGGCGGGGTCAGCGGATACGGCGTGGCTGGAGCAAGCT'
            'AGGGGGCAATGGAATTCCCGGTGTAGTGGTGGAATGCGTAGATATCGGGAGGAACACCAG'
            'TGGCGAAGGCGGTTGCCTGGAGCTTTGCTGACGCTGAGGCGCGAAAGCGTGGGGAGCGAT'
            'CCGGATTAGATACCCGGGTAGTCCACGCCGTAAACGATGCGGACTAGGTGTCGGGGGTAT'
            'CGACCCCCTCGGCGCCGCAGCTAACGCATTAAGTCCGCCGCCTGGGGACTACGGCCGCAA'
            'GGCTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCAGCGGAGCGTGTGGTTTAATT'
            'CGATGCAACGCGCAGAACCTTACCTCGGCTTGACATGCACCTGGTACCGAGGGGAAACCT'
            'GAGGGACCCGCAAGGGAGGGTGCACAGATGCTGCATGGCTGTCGTCAGCTCGTGCCGTGA'
            'GGTGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCTTAGTTGCCATCTCTAGGGA'
            'GACCGCCGGTCTAAACCGGAGGAAGGTGGGGATGACGTCAAGTCAGCACGGCTCTTACGT'
            'CGAGGGCTACACACACGCTACAATGGCCGGTACAACGGGCTGCAAAGGGGCGACCTGGAG'
            'CTAATCCCCCAAAGCCGGTCCCAGTTCGGATCGCAGCTGCAACTCCCCCTGCGTGAAGTC'
            'GGAGTTGCTAATATATCGGCTGGGTCAGCCCCCCGCAGGGATTTCGTTTCCCCGCCC'
        ),
        regions=[],
        rna_type='SO:0001000',
        url='https://www.arb-silva.de/browser/ssu/FN662328',
        seq_version='1',
        inference='Bacteria; Chloroflexi; TK10',
        common_name=None,
        species='uncultured eukaryote',
        lineage='Eukaryota; environmental samples; uncultured eukaryote',
        references=[
            pubs.reference('doi:10.1093/nar/gks1219'),
        ],
        description = 'uncultured eukaryote bacterial SSU rRNA',
    )
