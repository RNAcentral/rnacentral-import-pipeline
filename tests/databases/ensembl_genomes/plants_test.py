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


def parse(filename):
    with open(filename, 'r') as raw:
        return list(parser.parse(raw))


def entries_for(entries, accession):
    return [e for e in entries if e.accession == accession]


def entry_for(entries, accession):
    val = entries_for(entries, accession)
    assert len(val) == 1
    return val[0]


def has_entry_for(entries, accession):
    return bool(entries_for(entries, accession))


@pytest.fixture(scope='module')  # pylint: disable=no-member
def cress_2():
    return parse('data/ensembl_plants/Arabidopsis_thaliana.TAIR10.40.chromosome.2.dat')


@pytest.fixture(scope='module')  # pylint: disable=no-member
def oryza_9():
    return parse('data/ensembl_plants/Oryza_barthii.O.barthii_v1.41.chromosome.9.dat')


@pytest.fixture(scope='module')  # pylint: disable=no-member
def hordeum_pt():
    return parse('data/ensembl_plants/Hordeum_vulgare.IBSC_v2.41.chromosome.Pt.dat')


@pytest.fixture(scope='module')  # pylint: disable=no-member
def zea_7():
    return parse('data/ensembl_plants/Zea_mays.B73_RefGen_v4.41.chromosome.7.dat')


def test_can_parse_data(cress_2):
    val = attr.asdict(entry_for(cress_2, 'ENSEMBL_PLANTS:AT2G01010.1'))
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
                coordinate_system=dat.CoordinateSystem.one_based(),
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


def test_can_create_tair_entry(cress_2):
    val = attr.asdict(entry_for(cress_2, 'TAIR:AT2G01010.1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='AT2G01010.1',
        accession='TAIR:AT2G01010.1',
        ncbi_tax_id=3702,
        database='TAIR',
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
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='rRNA',
        url='',
        seq_version='1',
        note_data={},
        xref_data={
            'RefSeq': ['NR_139968.1'],
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


@pytest.mark.parameterize('accession,status', [
    ('TAIR:AT2G01010.1', True),
    ('TAIR:ENSRNA049757808-T1', False),
])
def generates_expected_inferred_entries(cress_2, accession, status):
    assert has_entry_for(cress_2, accession) == status


def test_can_get_with_odd_rna_type(cress_2):
    val = attr.asdict(entry_for(cress_2, 'ENSEMBL_PLANTS:AT2G03895.1'))
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
                coordinate_system=dat.CoordinateSystem.one_based(),
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
    val = attr.asdict(entry_for(cress_2, 'ENSEMBL_PLANTS:ENSRNA049492366-T1'))
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
                coordinate_system=dat.CoordinateSystem.one_based(),
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
    assert attr.asdict(entry_for(cress_2, 'ENSEMBL_PLANTS:AT2G03905.1')) == attr.asdict(dat.Entry(
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
                coordinate_system=dat.CoordinateSystem.one_based(),
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
    val = attr.asdict(entry_for(cress_2, 'ENSEMBL_PLANTS:ENSRNA049757815-T1'))
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
                coordinate_system=dat.CoordinateSystem.one_based(),
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


# OBART09G00240.1  transposable_elements
# ENSRNA049475598-T1  sense_intronic
@pytest.mark.skip()
def test_skips_transposable_elements(oryza_9):
    pass


def test_can_parse_rice_trna(oryza_9):
    val = attr.asdict(entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049456349-T1"))
    assert val == attr.asdict(dat.Entry(
        primary_id="ENSRNA049456349-T1",
        accession='ENSEMBL_PLANTS:ENSRNA049456349-T1',
        ncbi_tax_id=65489,
        database='ENSEMBL_PLANTS',
        sequence='TCCGTTGTAGTCTAGCTGGTTAGGATACTCGGCTCTCACCCGAGAGACCCGGGTTCGAGTCCCGGCAACGGAA',
        regions=[
            dat.SequenceRegion(
                chromosome='9',
                strand=1,
                exons=[dat.Exon(start=747032, stop=747104)],
                assembly_id='O.barthii_v1',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='tRNA',
        url='',
        seq_version='1',
        species='Oryza barthii',
        common_name='African wild rice',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; '
            'Oryzeae; Oryzinae; Oryza; Oryza barthii'
        ),
        gene="ENSRNA049456349",
        locus_tag="tRNA-Glu",
        description='Oryza barthii (African wild rice) tRNA-Glu for anticodon CUC',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_rice_snorna(oryza_9):
    val = attr.asdict(entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049475670-T1"))
    assert val == attr.asdict(dat.Entry(
        primary_id="ENSRNA049475670-T1",
        accession='ENSEMBL_PLANTS:ENSRNA049475670-T1',
        ncbi_tax_id=65489,
        database='ENSEMBL_PLANTS',
        sequence='AAAAAAGCAGGATGCTGTGTTCTCTATAAGCAGTGTCCTCGTAAATTTTAGGAACATGTTTCATCGTTATTGGGTGAACCGTTGGGCTATTCAATGTCCATTGGTTCAGTAAATGATGGCACATTT',
        regions=[
            dat.SequenceRegion(
                chromosome='9',
                strand=-1,
                exons=[dat.Exon(start=3344023, stop=3344148)],
                assembly_id='O.barthii_v1',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='snoRNA',
        url='',
        seq_version='1',
        species='Oryza barthii',
        common_name='African wild rice',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; '
            'Oryzeae; Oryzinae; Oryza; Oryza barthii'
        ),
        gene='ENSRNA049475670',
        locus_tag="snoR74",
        description='Oryza barthii (African wild rice) small nucleolar RNA snoR74',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_rice_pre_mirna(oryza_9):
    val = attr.asdict(entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049475651-T1"))
    assert val == attr.asdict(dat.Entry(
        primary_id='ENSRNA049475651-T1',
        accession='ENSEMBL_PLANTS:ENSRNA049475651-T1',
        ncbi_tax_id=65489,
        database='ENSEMBL_PLANTS',
        sequence='CCTCCCCGCCGGACCTCCCAGTGAGGAGGCTAGGGCCGCCAGGTCCGGTGATCCCATTCTCCTTGCCGGCGGATTCTGCGCCCTAGA',
        regions=[
            dat.SequenceRegion(
                chromosome='9',
                strand=-1,
                exons=[dat.Exon(start=3622031, stop=3622117)],
                assembly_id='O.barthii_v1',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='precursor_RNA',
        url='',
        seq_version='1',
        species='Oryza barthii',
        common_name='African wild rice',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; '
            'Oryzeae; Oryzinae; Oryza; Oryza barthii'
        ),
        gene='ENSRNA049475651',
        locus_tag='MIR1846',
        description='Oryza barthii (African wild rice) microRNA MIR1846',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_rice_u6(oryza_9):
    val = attr.asdict(entry_for(oryza_9, 'ENSEMBL_PLANTS:ENSRNA049475710-T1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='ENSRNA049475710-T1',
        accession='ENSEMBL_PLANTS:ENSRNA049475710-T1',
        ncbi_tax_id=65489,
        database='ENSEMBL_PLANTS',
        sequence='GTAGCTTATATACGCTGCTGTGCATAAAATTGAAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGAAGACGCACACAAATCGAGAAGTGGTCCAAATTTTT',
        regions=[
            dat.SequenceRegion(
                chromosome='9',
                strand=1,
                exons=[dat.Exon(start=6092721, stop=6092827)],
                assembly_id='O.barthii_v1',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type="snRNA",
        url='',
        seq_version='1',
        species='Oryza barthii',
        common_name='African wild rice',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; '
            'Oryzeae; Oryzinae; Oryza; Oryza barthii'
        ),
        gene='ENSRNA049475710',
        locus_tag="U6",
        description='Oryza barthii (African wild rice) U6 spliceosomal RNA',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_barley_antisense(hordeum_pt):
    val = attr.asdict(entry_for(hordeum_pt, 'ENSEMBL_PLANTS:ENSRNA049483195-T1'))
    assert val == attr.asdict(dat.Entry(
        primary_id='ENSRNA049483195-T1',
        accession='ENSEMBL_PLANTS:ENSRNA049483195-T1',
        ncbi_tax_id=112509,
        database='ENSEMBL_PLANTS',
        sequence='AATAACCAAATATAACACTGGGACTAAGGGTCAAATTGGTAATTTTTCTTACATCTCCCCCCCCAGGGGCCCAGGTATCATATACACCGCCAAAATAAAGAGCCTTGAGTACTAGAAGAAAAGCACCTAGACCTAACAAAATTAAGTGAATACCCAAAATTGTAGTCATTTTATTTCTATCTTTCCA',
        regions=[
            dat.SequenceRegion(
                chromosome='Pt',
                strand=-1,
                exons=[dat.Exon(start=10076, stop=10262)],
                assembly_id='IBSC_v2',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='antisense_RNA',
        url='',
        seq_version='1',
        species='Hordeum vulgare subsp. vulgare',
        common_name='two-rowed barley',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; Embryophyta; '
            'Tracheophyta; Spermatophyta; Magnoliophyta; Liliopsida; '
            'Poales; Poaceae; BOP clade; Pooideae; Triticodae; '
            'Triticeae; Hordeinae; Hordeum; Hordeum vulgare subsp. vulgare'
        ),
        gene='ENSRNA049483195',
        locus_tag='IsrR',
        description='Hordeum vulgare subsp. vulgare (two-rowed barley) antisense RNA which regulates isiA expression',
        references=[pubs.reference(29092050)],
    ))


def test_can_parse_zea_lincrna(zea_7):
    val = attr.asdict(entry_for(zea_7, "ENSEMBL_PLANTS:Zm00001d001070_T001"))
    assert val == attr.asdict(dat.Entry(
        primary_id="Zm00001d001070_T001",
        accession="ENSEMBL_PLANTS:Zm00001d001070_T001",
        ncbi_tax_id=4577,
        database='ENSEMBL_PLANTS',
        sequence=(
            'GTATGGAACACGGCCGACCCCAGGCATTGGCTTCTGGAGGTTGAAGATGGGCGCATGTCC'
            'GAGCGATCGGATGTGAATGGCTGTGGATAGTTGCGTGGTAGTGGTGGATGGCCAATCACT'
            'GGCGTAGCCATCGCCCTGGGTGCAGAACGTGGTCCGTATGGGGTCAGCTATGGCGCCGCC'
            'GCGCCGGACCCTGTTCACCTCCGTGGTTGCGGCCAGTGTGGGAAGATGGGCGAGCGCCGT'
            'TGGTATGGCCTGGAGCGGCTAGGATTAGGTGAGCACCTGGGTTGGGCGGGTTAAGTCCTG'
            'GGCGGTTAGAT'
        ),
        regions=[
            dat.SequenceRegion(
                chromosome='7',
                strand=-1,
                exons=[dat.Exon(start=359423, stop=359733)],
                assembly_id='B73_RefGen_v4',
                coordinate_system=dat.CoordinateSystem.one_based(),
            )
        ],
        rna_type='lncRNA',
        url='',
        seq_version='1',
        species='Zea mays',
        common_name='maize',
        lineage=(
            'Eukaryota; Viridiplantae; Streptophyta; '
            'Embryophyta; Tracheophyta; Spermatophyta; Magnoliophyta; '
            'Liliopsida; Poales; Poaceae; PACMAD clade; Panicoideae; '
            'Andropogonodae; Andropogoneae; Tripsacinae; Zea; Zea mays'
        ),
        gene="Zm00001d001070",
        description='Zea mays (maize) lncRNA Zm00001d001070',
        references=[pubs.reference(29092050)],
    ))


def test_does_not_generate_tair_for_others(zea_7):
    assert has_entry_for(zea_7, "TAIR:Zm00001d001070_T001") is False
    assert has_entry_for(zea_7, "ENSEMBL_PLANTS:Zm00001d001070_T001") is True
