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

import attr
import pytest

from rnacentral_pipeline.databases import data as dat
from rnacentral_pipeline.databases.helpers import publications as pubs

from . helpers import parse, entry_for


@pytest.fixture(scope='module')  # pylint: disable=no-member
def cress_2():
    return parse('data/ensembl_plants/Arabidopsis_thaliana.TAIR10.40.chromosome.2.dat')


def test_can_parse_data(cress_2):
    val = attr.asdict(entry_for(cress_2, 'AT2G01010.1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='AT2G01010.1',
        accession='ENSEMBL_PLANTS:AT2G01010.1',
        ncbi_tax_id=3702,
        database='ENSEMBL_PLANTS',
        sequence=(
            'TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGT'
            'AAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTG'
            'TTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAA'
            'CCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTT'
            'GCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCAT'
            'TCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGG'
            'GTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCA'
            'AGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAAT'
            'AACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAAC'
            'GAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGC'
            'GTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGG'
            'TCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTT'
            'AATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAA'
            'GCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTT'
            'GGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGT'
            'CAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGAT'
            'GTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTC'
            'TCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCA'
            'CCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTT'
            'AAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACA'
            'CGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGA'
            'TTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTC'
            'CGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCG'
            'GCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGT'
            'GATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACA'
            'CCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGC'
            'AATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTA'
            'CGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGT'
            'GTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAA'
            'CCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAG'
            'GATCATTG'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='2',
                strand=1,
                exons=[dat.Exon(start=3706, stop=5513)],
                assembly_id='TAIR10',
            )
        ],
        rna_type='rRNA',
        url='',
        seq_version='1',
        note_data={},
        xref_data={
            'RefSeq': ['NR_139968.1'],
            'TAIR': ['AT2G01010.1'],
        },
        species='Arabidopsis thaliana',
        common_name='thale-cress',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; '
            'Brassicales; Brassicaceae; Camelineae; Arabidopsis; '
            'Arabidopsis thaliana'
        ),
        gene='AT2G01010',
        description='Arabidopsis thaliana (thale-cress) rRNA',
        references=[pubs.reference(29092050)],
    ))


def test_can_get_with_odd_rna_type(cress_2):
    val = attr.asdict(entry_for(cress_2, 'AT2G03895.1'))
    print(val['sequence'])
    assert val == attr.asdict(dat.Entry(
        primary_id='AT2G03895.1',
        accession='ENSEMBL_PLANTS:AT2G03895.1',
        ncbi_tax_id=3702,
        database='ENSEMBL_PLANTS',
        sequence=(
            'GGTGGTCTCTGTTGGTGAATCGTCGTCATTGAGAGCTGACACCGGCCCAAAGCCTTTGCT'
            'CCGGCGTTGCGTGACGGAGTATCGGAGTCCAGCTTCCCTCCACGAATTGCAGAAAGTTAC'
            'AGCGTAAGGACAACGCTGCTTTGTAGGCGAACCCAAGTTGCGAGTGGTGAGGCGGAAATG'
            'GTGGATAAGAGCAGAACTAGTGCTTGTGCTGCTC'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='2',
                strand=1,
                exons=[dat.Exon(start=84873, stop=85086)],
                assembly_id='TAIR10',
            )
        ],
        rna_type='lncRNA',
        url='',
        seq_version='1',
        note_data={},
        xref_data={
            'RefSeq': ['NR_139974.1'],
        },
        species='Arabidopsis thaliana',
        common_name='thale-cress',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; '
            'Brassicales; Brassicaceae; Camelineae; Arabidopsis; '
            'Arabidopsis thaliana'
        ),
        gene='AT2G03895',
        description='Arabidopsis thaliana (thale-cress) lncRNA AT2G03895',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_a_trna(cress_2):
    val = attr.asdict(entry_for(cress_2, 'ENSRNA049492366-T1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='ENSRNA049492366-T1',
        accession='ENSEMBL_PLANTS:ENSRNA049492366-T1',
        ncbi_tax_id=3702,
        database='ENSEMBL_PLANTS',
        sequence=(
            'GCTGGAATAGCTCAGTTGGTTAGAGCGTGTGGCTGTTAACCACAAGGTCGGAGGTTCGAC'
            'CCCTCCTTCTAGCG'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='2',
                strand=1,
                exons=[dat.Exon(start=102065, stop=102138)],
                assembly_id='TAIR10',
            )
        ],
        rna_type='tRNA',
        url='',
        seq_version='1',
        species='Arabidopsis thaliana',
        common_name='thale-cress',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; '
            'Brassicales; Brassicaceae; Camelineae; Arabidopsis; '
            'Arabidopsis thaliana'
        ),
        gene='ENSRNA049492366',
        locus_tag='tRNA-Asn',
        description='Arabidopsis thaliana (thale-cress) tRNA-Asn for anticodon GUU',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_gene_with_minimal_metadata(cress_2):
    assert attr.asdict(entry_for(cress_2, 'AT2G03905.1')) == attr.asdict(dat.Entry(
        primary_id='AT2G03905.1',
        accession='ENSEMBL_PLANTS:AT2G03905.1',
        ncbi_tax_id=3702,
        database='ENSEMBL_PLANTS',
        sequence=(
            'CGCCGTTAGTCCGTGAGGAGAAAATAGGCCCACTCTGGCACACTCTCTCTGGGTTTAGGT'
            'TTAGGTTTTTTTGGGGCTCTCTATCCTAAGAAACTAGGAGACATCACACTTCACCAAGTC'
            'TACTTATCGACAATTTTATCGTATCACCATAACGACAATAAGGGCCGGACTAATGTTTGT'
            'ACACATGTCCTCTCCTTTTACCCTT'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='2',
                strand=-1,
                exons=[dat.Exon(start=122882, stop=123086)],
                assembly_id='TAIR10',
            )
        ],
        rna_type='lncRNA',
        url='',
        seq_version='1',
        species='Arabidopsis thaliana',
        common_name='thale-cress',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; '
            'Brassicales; Brassicaceae; Camelineae; Arabidopsis; '
            'Arabidopsis thaliana'
        ),
        gene='AT2G03905',
        description='Arabidopsis thaliana (thale-cress) lncRNA AT2G03905',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_premirna(cress_2):
    val = attr.asdict(entry_for(cress_2, 'ENSRNA049757815-T1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='ENSRNA049757815-T1',
        accession='ENSEMBL_PLANTS:ENSRNA049757815-T1',
        ncbi_tax_id=3702,
        database='ENSEMBL_PLANTS',
        sequence=(
            'CGTAAAGCAGGTGATTCACCAATTTAGGTTTACATCCACAGTGTGGAAGACACTGAAGGA'
            'CCTAAACTAACAAAGGTAAACGGCTCAGTGTGCGGGGTATTACACTCGGTTTAATGTCTG'
            'AATGCGATAATCCGCACGATGATCTCTTTATCTTTGTTTGTTTAGGTCCCTTAGTTTCTT'
            'CTATACCGTGAATCCAATCCTGTATTGGATGAGCTGGTTTATGACC'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='2',
                strand=-1,
                exons=[dat.Exon(start=771337, stop=771562)],
                assembly_id='TAIR10',
            )
        ],
        rna_type='precursor_RNA',
        url='',
        seq_version='1',
        species='Arabidopsis thaliana',
        common_name='thale-cress',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; '
            'Brassicales; Brassicaceae; Camelineae; Arabidopsis; '
            'Arabidopsis thaliana'
        ),
        gene='ENSRNA049757815',
        locus_tag='ath-MIR840',
        description='Arabidopsis thaliana (thale-cress) precursor RNA ath-MIR840',
        references=[pubs.reference(29092050)],
    ))
