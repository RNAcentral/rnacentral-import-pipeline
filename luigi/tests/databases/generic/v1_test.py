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

import json

import attr
import pytest

from databases.generic import v1


@pytest.mark.parametrize('filename,taxids', [
    ('data/json-schema/v020/flybase.json', [7227, 7227, 7227, 7227, 7227]),
    ('data/json-schema/v020/lincipedia.json', [9606]),
])
def test_can_extract_taxid(filename, taxids):
    with open(filename, 'r') as raw:
        data = json.load(raw)['data']
        assert [v1.taxid(e) for e in data] == taxids


@pytest.mark.parametrize('filename,xrefs', [
    ('data/json-schema/v020/lincipedia.json', [{"NONCODE": ["NONHSAT050743"]}]),
])
def test_can_generate_xref_data(filename, xrefs):
    with open(filename, 'r') as raw:
        data = json.load(raw)['data']
        assert [v1.xrefs(e) for e in data] == xrefs


@pytest.mark.skip()
def test_can_extract_anticodon():
    pass


@pytest.mark.parametrize('filename,count', [
    ('data/json-schema/v020/flybase.json', 5),
    ('data/json-schema/v020/lincipedia.json', 1),
])
def test_can_parse_all_data(filename, count):
    with open(filename, 'r') as raw:
        data = json.load(raw)
        assert len(list(v1.parse(data))) == count


def test_can_correctly_parse_data():
    with open('data/json-schema/v020/flybase.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert attr.asdict(data[0]) == {
        'primary_id': 'FBtr0346876',
        'accession': 'FLYBASE:FBtr0346876',
        'ncbi_tax_id': 7227,
        'database': 'FLYBASE',
        'sequence': 'TTATATACAACCTCAACTCATATGGGACTACCCCCTGAATTTAAGCATATTAATTAGGGGAGGAAAAGAAACTAACAAGGATTTTCTTAGTAGCGGCGAGCGAAAAGAAAACAGTTCAGCACTAAGTCACTTTGTCTATATGGCAAATGTGAGATGCAGTGTATGGAGCGTCAATATTCTAGTATGAGAAATTAACGATTTAAGTCCTTCTTAAATGAGGCCATTTACCCATAGAGGGTGCCAGGCCCGTATAACGTTAATGATTACTAGATGATGTTTCCAAAGAGTCGTGTTGCTTGATAGTGCAGCACTAAGTGGGTGGTAAACTCCATCTAAAACTAAATATAACCATGAGACCGATAGTAAACAAGTACCGTGAGGGAAAGTTGAAAAGAACTCTGAATAGAGAGTTAAACAGTACGTGAAACTGCTTAGAGGTTAAGCCCGATGAACCTGAATATCCGTTATGGAAAATTCATCATTAAAATTGTAATATTTAAATAATATTATGAGAATAGTGTGCATTTTTTCCATATAAGGACATTGTAATCTATTAGCATATACCAAATTTATCATAAAATATAACTTATAGTTTATTCCAATTAAATTGCTTGCATTTTAACACAGAATAAATGTTATTAATTTGATAAAGTGCTGATAGATTTATATGATTACAGTGCGTTAATTTTTCGGAATTATATAATGGCATAATTATCATTGATTTTTGTGTTTATTATATGCACTTGTATGATTAACAATGCGAAAGATTCAGGATACCTTCGGGACCCGTCTTGAAACACGGACCAAGGAGTCTAACATATGTGCAAGTTATTGGGATATAAACCTAATAGCGTAATTAACTTGACTAATAATGGGATTAGTTTTTTAGCTATTTATAGCTGCTAATTAACACAATCCCGGGGCGTTCTATATAGTTATGTATAATGTATATTTATATTATTTATGCCTCTAACTGGAACGTACCTTGAGCATATATGCTGTGACCCGAAAGATGGTGAACTATACTTGATCAGGTTGAAGTCAGGGGAAACCCTGATGGAAGACCGAAACAGTTCTGACGTGCAAATCGATTGTCAGAATTGAGTATAGGGGCGAAAGACCAATCGAACCATCTAGTAGCTGGTTCCTTCCGAAGTTTCCCTCAGGATAGCTGGTGCATTTTAATATTATATAAAATAATCTTATCTGGTAAAGCGAATGATTAGAGGCCTTAGGGTCGAAACGATCTTAACCTATTCTCAAACTTTAAATGGGTAAGAACCTTAACTTTCTTGATATGAAGTTCAAGGTTATGATATAATGTGCCCAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGTGGGATGAACCAAACGTAATGTTACGGTGCCCAAATTAACAACTCATGCAGATACCATGAAAGGCGTTGGTTGCTTAAAACAGCAGGACGGTGATCATGGAAGTCGAAATCCGCTAAGGAGTGTGTAACAACTCACCTGCCGAAGCAACTAGCCCTTAAAATGGATGGCGCTTAAGTTGTATACCTATACATTACCGCTAAAGTAGATGATTTATATTACTTGTGATATAAATTTTGATCTTGGTGGTAGTAGCAAATAATCGAATGAGACCTTGGAGGACTGAAGTGGAGAAGGGTTTCGTGTGAACAGTGGTTGATCACGAGTTAGTCGGTCCTAAGTTCAAGGCGAAAGCCGAAAATTTTCAAGTAAAACAAAAATGCCTAACTATATAAACAAAGCGAATTATAATACACTTGAATAATTTTGAACGAAAGGGAATACGGTTCCAATTCCGTAACCTGTTGAGTATCCGTTTGTTTTTCTGTTTTATAGCCGTACTACCATGGAAGTCTTTCGCAGAGAGATATGGTAGATGGGCTAGAAGAGCATGACATATACTGTTGTGTCGATATTTTCTCCTCGGACCTTGAAAATTTATGGTGGGGACACGCAAACTTCTCAACAGGCCGTACCAATATCCGCAGCTGGTCTCCAAGGTGAAGAGTCTCTAGTCGATAGAATAATGTAGGTAAGGGAAGTCGGCAAATTAGATCCGTAACTTCGGGATAAGGATTGGCTCTGAAGATTGAGATAGTCGGGCTTGATTGGGAAACAATAACATGGTTTATGTGCTCGTTTCTTGCTAAAATTTTTAAGAATACTATTTGGGTTAAACCAATTAGTTCTTATTAATTATAACGATTATCAATTAACAATCAATTCAGAACTGGCACGGACTTGGGGAATCCGACTGTCTAATTAAAACAAAGCATTGTGATGGCCCTAGCGGGTGTTGACACAATGTGATTTCTGCCCAGTGCTCTGAATGTCAAAGTGAAGAAATTCAAGTAAGCGCGGGTCAACGGCGGGAGTAACTATGACTCTCTTAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATTAACGAGATTCCTACT',
        'exons': [{
            'chromosome': 'rDNA',
            'complement': False,
            'primary_end': 49485,
            'primary_start': 46771,
        }],
        'rna_type': 'SO:0000252',
        'url': 'http://flybase.org/reports/FBtr0346876.html',
        'seq_version': '1',
        'note_data': {},
        'xref_data': {
            'REFSEQ': ['NR_133553.1'],
        },
        'chromosome': None,
        'species': 'Drosophila melanogaster',
        'common_name': None,
        'lineage': (
            'Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; '
            'Insecta; Pterygota; Neoptera; Holometabola; Diptera; Brachycera; '
            'Muscomorpha; Ephydroidea; Drosophilidae; Drosophila; Sophophora; '
            'Drosophila melanogaster'
        ),
        'gene': 'Dmel_CR45837',
        'locus_tag': None,
        'optional_id': None,
        'product': '28S ribosomal RNA',
        'parent_accession': 'CP007120',
        'ordinal': None,
        'non_coding_id': None,
        'project': None,
        'keywords': None,
        'division': None,
        'organelle': None,
        'allele': None,
        'anticodon': None,
        'experiment': None,
        'function': None,
        'inference': None,
        'map': None,
        'old_locus_tag': None,
        'operon': None,
        'standard_name': None,
        'description': 'Drosophila melanogaster (fruit fly) 28S ribosomal RNA:CR45837',
        'mol_type': None,
        'is_composite': None,
        'pseudogene': None,
        'gene_synonyms': [
            "CR45837",
        ],
        'references': [
        ],
        'secondary_structure': {'dot_bracket': ''},
    }
