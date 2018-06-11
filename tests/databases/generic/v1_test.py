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

from databases import data as dat
from databases.helpers import publications as pub

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
        'sequence': "TTATATACAACCTCAACTCATATGGGACTACCCCCTGAATTTAAGCATATTAATTAGGGGAGGAAAAGAAACTAACAAGGATTTTCTTAGTAGCGGCGAGCGAAAAGAAAACAGTTCAGCACTAAGTCACTTTGTCTATATGGCAAATGTGAGATGCAGTGTATGGAGCGTCAATATTCTAGTATGAGAAATTAACGATTTAAGTCCTTCTTAAATGAGGCCATTTACCCATAGAGGGTGCCAGGCCCGTATAACGTTAATGATTACTAGATGATGTTTCCAAAGAGTCGTGTTGCTTGATAGTGCAGCACTAAGTGGGTGGTAAACTCCATCTAAAACTAAATATAACCATGAGACCGATAGTAAACAAGTACCGTGAGGGAAAGTTGAAAAGAACTCTGAATAGAGAGTTAAACAGTACGTGAAACTGCTTAGAGGTTAAGCCCGATGAACCTGAATATCCGTTATGGAAAATTCATCATTAAAATTGTAATATTTAAATAATATTATGAGAATAGTGTGCATTTTTTCCATATAAGGACATTGTAATCTATTAGCATATACCAAATTTATCATAAAATATAACTTATAGTTTATTCCAATTAAATTGCTTGCATTTTAACACAGAATAAATGTTATTAATTTGATAAAGTGCTGATAGATTTATATGATTACAGTGCGTTAATTTTTCGGAATTATATAATGGCATAATTATCATTGATTTTTGTGTTTATTATATGCACTTGTATGATTAACAATGCGAAAGATTCAGGATACCTTCGGGACCCGTCTTGAAACACGGACCAAGGAGTCTAACATATGTGCAAGTTATTGGGATATAAACCTAATAGCGTAATTAACTTGACTAATAATGGGATTAGTTTTTTAGCTATTTATAGCTGCTAATTAACACAATCCCGGGGCGTTCTATATAGTTATGTATAATGTATATTTATATTATTTATGCCTCTAACTGGAACGTACCTTGAGCATATATGCTGTGACCCGAAAGATGGTGAACTATACTTGATCAGGTTGAAGTCAGGGGAAACCCTGATGGAAGACCGAAACAGTTCTGACGTGCAAATCGATTGTCAGAATTGAGTATAGGGGCGAAAGACCAATCGAACCATCTAGTAGCTGGTTCCTTCCGAAGTTTCCCTCAGGATAGCTGGTGCATTTTAATATTATATAAAATAATCTTATCTGGTAAAGCGAATGATTAGAGGCCTTAGGGTCGAAACGATCTTAACCTATTCTCAAACTTTAAATGGGTAAGAACCTTAACTTTCTTGATATGAAGTTCAAGGTTATGATATAATGTGCCCAGTGGGCCACTTTTGGTAAGCAGAACTGGCGCTGTGGGATGAACCAAACGTAATGTTACGGTGCCCAAATTAACAACTCATGCAGATACCATGAAAGGCGTTGGTTGCTTAAAACAGCAGGACGGTGATCATGGAAGTCGAAATCCGCTAAGGAGTGTGTAACAACTCACCTGCCGAAGCAACTAGCCCTTAAAATGGATGGCGCTTAAGTTGTATACCTATACATTACCGCTAAAGTAGATGATTTATATTACTTGTGATATAAATTTTGAAACTTTAGTGAGTAGGAAGGTACAATGGTATGCGTAGAAGTGTTTGGCGTAAGCCTGCATGGAGCTGCCATTGGTACAGATCTTGGTGGTAGTAGCAAATAATCGAATGAGACCTTGGAGGACTGAAGTGGAGAAGGGTTTCGTGTGAACAGTGGTTGATCACGAGTTAGTCGGTCCTAAGTTCAAGGCGAAAGCCGAAAATTTTCAAGTAAAACAAAAATGCCTAACTATATAAACAAAGCGAATTATAATACACTTGAATAATTTTGAACGAAAGGGAATACGGTTCCAATTCCGTAACCTGTTGAGTATCCGTTTGTTATTAAATATGGGCCTCGTGCTCATCCTGGCAACAGGAACGACCATAAAGAAGCCGTCGAGAGATATCGGAAGAGTTTTCTTTTCTGTTTTATAGCCGTACTACCATGGAAGTCTTTCGCAGAGAGATATGGTAGATGGGCTAGAAGAGCATGACATATACTGTTGTGTCGATATTTTCTCCTCGGACCTTGAAAATTTATGGTGGGGACACGCAAACTTCTCAACAGGCCGTACCAATATCCGCAGCTGGTCTCCAAGGTGAAGAGTCTCTAGTCGATAGAATAATGTAGGTAAGGGAAGTCGGCAAATTAGATCCGTAACTTCGGGATAAGGATTGGCTCTGAAGATTGAGATAGTCGGGCTTGATTGGGAAACAATAACATGGTTTATGTGCTCGTTCTGGGTAAATAGAGTTTCTAGCATTTATGTTAGTTACTTGTTCCCCGGATAGTTTAGTTACGTAGCCAATTGTGGAACTTTCTTGCTAAAATTTTTAAGAATACTATTTGGGTTAAACCAATTAGTTCTTATTAATTATAACGATTATCAATTAACAATCAATTCAGAACTGGCACGGACTTGGGGAATCCGACTGTCTAATTAAAACAAAGCATTGTGATGGCCCTAGCGGGTGTTGACACAATGTGATTTCTGCCCAGTGCTCTGAATGTCAAAGTGAAGAAATTCAAGTAAGCGCGGGTCAACGGCGGGAGTAACTATGACTCTCTTAAGGTAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATTAACGAGATTCCTACT",
        'exons': [{
            'chromosome_name': 'rDNA',
            'primary_start': 46772,
            'primary_end': 49485,
            "assembly_id": "R6",
            'complement': False,
        }],
        'rna_type': 'rRNA',
        'url': 'http://flybase.org/reports/FBtr0346876.html',
        'seq_version': '1',
        'note_data': {},
        'xref_data': {
            'REFSEQ': ['NR_133553.1'],
        },
        'chromosome': None,
        'species': 'Drosophila melanogaster',
        'common_name': 'fruit fly',
        'lineage': (
            'Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; '
            'Insecta; Pterygota; Neoptera; Holometabola; Diptera; Brachycera; '
            'Muscomorpha; Ephydroidea; Drosophilidae; Drosophila; Sophophora; '
            'Drosophila melanogaster'
        ),
        'gene': 'FBgn0267497',
        'locus_tag': 'Dmel_CR45837',
        'optional_id': None,
        'product': None,
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
        'location_start': None,
        'location_end': None,
        'gene_synonyms': [
            "CR45837",
        ],
        'references': [
        ],
        'secondary_structure': {'dot_bracket': ''},
    }


def test_can_correctly_parse_lncipedia_data():
    with open('data/json-schema/v020/lncipedia-5.0.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 2
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='lnc-CLEC18B-3:5',
        accession='LNCIPEDIA:lnc-CLEC18B-3:5',
        ncbi_tax_id=9606,
        database='LNCIPEDIA',
        sequence=(
            "GTAGATCATCATCATAACAGCTCCCAGTGAATCAATGCCTCCCTGCATCCACACCCCTT"
            "TGTAACATGATATTGTTGCTCTTCCCATCAAGAAATGGTCTCTTTTGGCCGGGCACAGTA"
            "GCTCACGTCATCCCAGCACTTTGGGAGGCTGAGGCAGGCAGATGGCTTGAGTGTAACAAT"
            "ATGTTGAATTTGCCATGGGCCTTTAAAGCTTCAATGTTGTGAAGAGCTCTGCATGAAATT"
            "TTAAAGAGACTGGACCTTTCATCTGCACAACAGCAGGGCACCTCGCTATAGGGACAAGAA"
            "AGGAAGAGAGAGAGAGAACATTTCTGAAGTAATAGTGAAAAATAACAGCAGAAGCAATTA"
            "TTTCATCAAAGATTGCAGGGAGAGGCTCCTCCGTGCTCCTGAGAGGCCGAACACAGGGTC"
            "GCCAGCACAGCTTACTGCTCGGTGTCTCCTGAGCCACAGAGGAAGACGTGGCAGGAGCAC"
            "CTGGTGCTAATATATATTCATGTCTATGGCAATGCCGACCATCTGGCTGGTCTGAACCAG"
            "GATAAAAGTGAAGAATTCCTCTGTGAAGACCCAGCTCTTTCTTTGGCTCCTTTTTTGAAG"
            "CCATCTTTGCTCTGCTCTCCTCTGCTGCCCAGAAAGTTCCAGAGTGAAGCTCAGCTCTAG"
            "ATGAACAAAAACTGGTTGAGTCCAGAGATGCCTGAGTTGGAGATGAACCTTGCAAACTTT"
            "CCTCATTACCATACTAAAAACCCCACCCAGGAAGGAGCTTATCTGCCATTTCCTACACAT"
            "GTGACATATGGAGAAGCATGATCAGCTACTTCACAGTCTCTGCCTTTACTCTGCCTCCGC"
            "ATACAATGGCTCAGCCAACTAGCCTAACGAAAGCTGTTTTCACCATTGTTTGGGAGGTAC"
            "TGCTTTGGGAAACTGCCCCAGCTGTCCTCCTTACTTGTTGTAGGTAATAAAATCCCTTTG"
            "TTAAATC"
        ),
        exons=[
            dat.Exon(chromosome_name='16', primary_start=74226291, primary_end=74226625, assembly_id='GRCh37', complement=True),
            dat.Exon(chromosome_name='16', primary_start=74239804, primary_end=74240064, assembly_id='GRCh37', complement=True),
            dat.Exon(chromosome_name='16', primary_start=74244205, primary_end=74244404, assembly_id='GRCh37', complement=True),
            dat.Exon(chromosome_name='16', primary_start=74249251, primary_end=74249420, assembly_id='GRCh37', complement=True),
        ],
        rna_type='SO:0001877',
        url='https://lncipedia.org/db/transcript/lnc-CLEC18B-3:5',
        seq_version='1',

        xref_data={'NONCODE': ['NONHSAT143655']},

        gene='lnc-CLEC18B-3',
        gene_synonyms=[
            "ENSG00000249447",
            "XLOC_012007",
            "linc-ZFHX3-2",
            "ENSG00000261404.1",
            "AC009120.4",
            "OTTHUMG00000176255.2",
            "ENSG00000261404.5",
            "ENSG00000261404.6",
            "AC138627.1",
            "LOC101928035"
        ],
        # product='long non-coding RNA lnc-CLEC18B-3:5',

        description='Homo sapiens (human) non-protein coding lnc-CLEC18B-3:5',
        species='Homo sapiens',
        common_name='human',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Primates; Haplorrhini; Catarrhini; '
            'Hominidae; Homo; Homo sapiens'
        ),
    ))

def test_can_correctly_parse_mirbase_data():
    with open('data/json-schema/v020/missing-mirbase.json', 'r') as raw:
        data = json.load(raw)
        data = list(v1.parse(data))

    assert len(data) == 1
    assert attr.asdict(data[0]) == attr.asdict(dat.Entry(
        primary_id='MI0000612',
        accession='MIRBASE:MI0000612',
        ncbi_tax_id=10116,
        database='MIRBASE',
        sequence=(
            "TCTTTTGGGCGGGGGTCAAGAGCAATAACGAAAAATGTTTGTTTTTCGTAAACCGTTTTT"
            "CATTATTGCTCCTGACCTCCTCTCATTTGTTATAGCCA"
        ),
        exons=[],
        rna_type='SO:0001244',
        url='http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000612',
        seq_version='1',
        xref_data={
            'EntrezGene': ['Mir335'],
        },
        description='Rattus norvegicus miR-335 stem-loop',
        species='Rattus norvegicus',
        common_name='Norway rat',
        lineage=(
            'Eukaryota; Metazoa; Chordata; Craniata; '
            'Vertebrata; Euteleostomi; Mammalia; Eutheria; '
            'Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; '
            'Muridae; Murinae; Rattus; Rattus norvegicus'
        ),
        references=[
            pub.reference(17604727),
            pub.reference(14691248),
        ],
    ))
