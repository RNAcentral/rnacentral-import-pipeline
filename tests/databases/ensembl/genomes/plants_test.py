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
from rnacentral_pipeline.databases.ensembl import plants
from rnacentral_pipeline.databases.helpers import publications as pubs

from . import helpers


@pytest.fixture(scope="module")  # pylint: disable=no-member
def cress_2():
    return helpers.parse(
        plants.parse,
        "test-data/ensembl_plants/Arabidopsis_thaliana.TAIR10.chromosome.2.dat",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def oryza_9():
    return helpers.parse(
        plants.parse,
        "test-data/ensembl_plants/Oryza_barthii.O.barthii_v1.chromosome.9.dat",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def hordeum_pt():
    return helpers.parse(
        plants.parse,
        "test-data/ensembl_plants/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.nonchromosomal.dat",
    )


@pytest.fixture(scope="module")  # pylint: disable=no-member
def zea_7():
    return helpers.parse(
        plants.parse,
        "test-data/ensembl_plants/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.chromosome.7.dat",
    )


def test_can_parse_data(cress_2):
    val = attr.asdict(helpers.entry_for(cress_2, "ENSEMBL_PLANTS:AT2G01010.1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="AT2G01010.1",
            accession="ENSEMBL_PLANTS:AT2G01010.1",
            ncbi_tax_id=3702,
            database="ENSEMBL_PLANTS",
            sequence=(
                "TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGT"
                "AAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTG"
                "TTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAA"
                "CCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTT"
                "GCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCAT"
                "TCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGG"
                "GTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCA"
                "AGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAAT"
                "AACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAAC"
                "GAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGC"
                "GTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGG"
                "TCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTT"
                "AATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAA"
                "GCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTT"
                "GGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGT"
                "CAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGAT"
                "GTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTC"
                "TCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCA"
                "CCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTT"
                "AAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACA"
                "CGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGA"
                "TTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTC"
                "CGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCG"
                "GCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGT"
                "GATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACA"
                "CCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGC"
                "AATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTA"
                "CGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGT"
                "GTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAA"
                "CCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAG"
                "GATCATTG"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=1,
                    exons=[dat.Exon(start=31, stop=1838)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={
                "EntrezGene_trans_name": ["AT2G01010-201"],
                "RefSeq": ["NR_139968.1"],
                "RefSeq_ncRNA": ["NR_139968"],
                "TAIR": ["AT2G01010.1"],
                "GO": ["GO:0003735", "GO:0005840"],
                "RNAcentral": ["URS000008172F"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G01010",
            locus_tag="AT2G01010",
            description="Arabidopsis thaliana (thale-cress) 18S ribosomal RNA",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_create_tair_entry(cress_2):
    val = attr.asdict(helpers.entry_for(cress_2, "TAIR:AT2G01010.1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="AT2G01010.1",
            accession="TAIR:AT2G01010.1",
            ncbi_tax_id=3702,
            database="TAIR",
            sequence=(
                "TACCTGGTTGATCCTGCCAGTAGTCATATGCTTGTCTCAAAGATTAAGCCATGCATGTGT"
                "AAGTATGAACGAATTCAGACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTG"
                "TTTGATGGTAACTACTACTCGGATAACCGTAGTAATTCTAGAGCTAATACGTGCAACAAA"
                "CCCCGACTTATGGAAGGGACGCATTTATTAGATAAAAGGTCGACGCGGGCTCTGCCCGTT"
                "GCTCTGATGATTCATGATAACTCGACGGATCGCATGGCCTCTGTGCTGGCGACGCATCAT"
                "TCAAATTTCTGCCCTATCAACTTTCGATGGTAGGATAGTGGCCTACCATGGTGGTAACGG"
                "GTGACGGAGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACATCCA"
                "AGGAAGGCAGCAGGCGCGCAAATTACCCAATCCTGACACGGGGAGGTAGTGACAATAAAT"
                "AACAATACTGGGCTCTTTCGAGTCTGGTAATTGGAATGAGTACAATCTAAATCCCTTAAC"
                "GAGGATCCATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCCAGCTCCAATAGC"
                "GTATATTTAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAACCTTGGGATGGGTCGGCCGG"
                "TCCGCCTTTGGTGTGCATTGGTCGGCTTGTCCCTTCGGTCGGCGATACGCTCCTGGTCTT"
                "AATTGGCCGGGTCGTGCCTCCGGCGCTGTTACTTTGAAGAAATTAGAGTGCTCAAAGCAA"
                "GCCTACGCTCTGGATACATTAGCATGGGATAACATCATAGGATTTCGATCCTATTGTGTT"
                "GGCCTTCGGGATCGGAGTAATGATTAACAGGGACAGTCGGGGGCATTCGTATTTCATAGT"
                "CAGAGGTGAAATTCTTGGATTTATGAAAGACGAACAACTGCGAAAGCATTTGCCAAGGAT"
                "GTTTTCATTAATCAAGAACGAAAGTTGGGGGCTCGAAGACGATCAGATACCGTCCTAGTC"
                "TCAACCATAAACGATGCCGACCAGGGATCAGCGGATGTTGCTTATAGGACTCCGCTGGCA"
                "CCTTATGAGAAATCAAAGTTTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTT"
                "AAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACA"
                "CGGGGAAACTTACCAGGTCCAGACATAGTAAGGATTGACAGACTGAGAGCTCTTTCTTGA"
                "TTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTC"
                "CGTTAATGAACGAGACCTCAGCCTGCTAACTAGCTACGTGGAGGCATCCCTTCACGGCCG"
                "GCTTCTTAGAGGGACTATGGCCGTTTAGGCCAAGGAAGTTTGAGGCAATAACAGGTCTGT"
                "GATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGTATTCAACGAGTTCACA"
                "CCTTGGCCGACAGGCCCGGGTAATCTTTGAAATTTCATCGTGATGGGGATAGATCATTGC"
                "AATTGTTGGTCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAGCTCGCGTTGACTA"
                "CGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAAGT"
                "GTTCGGATCGCGGCGACGTGGGTGGTTCGCCGCCCGCGACGTCGCGAGAAGTCCACTAAA"
                "CCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAG"
                "GATCATTG"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=1,
                    exons=[dat.Exon(start=31, stop=1838)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={
                "EntrezGene_trans_name": ["AT2G01010-201"],
                "RefSeq": ["NR_139968.1"],
                "RefSeq_ncRNA": ["NR_139968"],
                "GO": ["GO:0003735", "GO:0005840"],
                "RNAcentral": ["URS000008172F"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G01010",
            locus_tag="AT2G01010",
            description="Arabidopsis thaliana (thale-cress) 18S ribosomal RNA",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_get_with_odd_rna_type(cress_2):
    val = attr.asdict(helpers.entry_for(cress_2, "ENSEMBL_PLANTS:AT2G03895.1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="AT2G03895.1",
            accession="ENSEMBL_PLANTS:AT2G03895.1",
            ncbi_tax_id=3702,
            database="ENSEMBL_PLANTS",
            sequence=(
                "GGTGGTCTCTGTTGGTGAATCGTCGTCATTGAGAGCTGACACCGGCCCAAAGCCTTTGCT"
                "CCGGCGTTGCGTGACGGAGTATCGGAGTCCAGCTTCCCTCCACGAATTGCAGAAAGTTAC"
                "AGCGTAAGGACAACGCTGCTTTGTAGGCGAACCCAAGTTGCGAGTGGTGAGGCGGAAATG"
                "GTGGATAAGAGCAGAACTAGTGCTTGTGCTGCTC"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=1,
                    exons=[dat.Exon(start=1899, stop=2112)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={
                "EntrezGene_trans_name": ["AT2G03895-201"],
                "RefSeq": ["NR_139974.1"],
                "RefSeq_ncRNA": ["NR_139974"],
                "RNAcentral": ["URS0000A76CC2"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G03895",
            locus_tag="AT2G03895",
            description="Arabidopsis thaliana (thale-cress) ncRNA",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_a_trna(cress_2):
    # ENSRNA049492366 (the original gene here) doesn't exist in the current
    # annotation - re-pointed at a current tRNA gene.
    val = attr.asdict(helpers.entry_for(cress_2, "ENSEMBL_PLANTS:AT2G07752.1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="AT2G07752.1",
            accession="ENSEMBL_PLANTS:AT2G07752.1",
            ncbi_tax_id=3702,
            database="ENSEMBL_PLANTS",
            sequence=(
                "GTCCCTTTCGTCCAGTGGTTAGGACATCGTCTTTTCATGTCGAAGACACGGGTTCGATT"
                "CCCGTAAGGGATA"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=-1,
                    exons=[dat.Exon(start=2438, stop=2509)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "TAIR": ["AT2G07752.1"],
                "GO": ["GO:0006412", "GO:0030533"],
                "RNAcentral": ["URS00002FBDA1"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G07752",
            description="Arabidopsis thaliana (thale-cress) misc RNA AT2G07752",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_gene_with_minimal_metadata(cress_2):
    assert attr.asdict(
        helpers.entry_for(cress_2, "ENSEMBL_PLANTS:AT2G03905.1")
    ) == attr.asdict(
        dat.Entry(
            primary_id="AT2G03905.1",
            accession="ENSEMBL_PLANTS:AT2G03905.1",
            ncbi_tax_id=3702,
            database="ENSEMBL_PLANTS",
            sequence=(
                "CGCCGTTAGTCCGTGAGGAGAAAATAGGCCCACTCTGGCACACTCTCTCTGGGTTTAGGT"
                "TTAGGTTTTTTTGGGGCTCTCTATCCTAAGAAACTAGGAGACATCACACTTCACCAAGTC"
                "TACTTATCGACAATTTTATCGTATCACCATAACGACAATAAGGGCCGGACTAATGTTTGT"
                "ACACATGTCCTCTCCTTTTACCCTT"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=-1,
                    exons=[dat.Exon(start=2173, stop=2377)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={
                "EntrezGene_trans_name": ["AT2G03905-201"],
                "RefSeq_ncRNA": ["NR_139975"],
                "RNAcentral": ["URS0000A7703F"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G03905",
            locus_tag="AT2G03905",
            description="Arabidopsis thaliana (thale-cress) ncRNA",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_premirna(cress_2):
    # ENSRNA049757815 (the original gene here) doesn't exist in the current
    # annotation - re-pointed at a current miRNA-precursor gene.
    val = attr.asdict(helpers.entry_for(cress_2, "ENSEMBL_PLANTS:at2g08015"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="at2g08015",
            accession="ENSEMBL_PLANTS:at2g08015",
            ncbi_tax_id=3702,
            database="ENSEMBL_PLANTS",
            sequence=(
                "AGTAGTCCACTGTGGTCTAAGGCTGAAATGCGGTGACCAGACTAGAAAACATGTCAATG"
                "GTGGTCCTTGTGAAAGGGTTGACAAAGAAGTCTAGGGAGACGGTGACACGTCGGATCTC"
                "TGGCCTATGTGAGTGCTAAAACAATTATGGTATTGGTTGCTTAAAGGTGTCATATTGCC"
                "ACTACGACTTTGAGAATGGTGAATGGTTGTGGCAATATGACTTTAAGCAGCCAATAAT"
                "ATAATCTTTTTAGCACTGATGGTGGCTAGAGAACCGGCATGCCACTGTCTCCGTAGAT"
                "TTTTTTTGTCAACCCTTTTGCAAGTCCCACCATTGATCAATTTTCCGAGCGGGTCACC"
                "GTATTTTAGCCTCAGACCACGGTGGACTACT"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="2",
                    strand=1,
                    exons=[dat.Exon(start=2570, stop=2951)],
                    assembly_id="TAIR10",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "RefSeq": ["NR_140265.1"],
                "RefSeq_ncRNA": ["NR_140265"],
                "miRBase_trans_name": ["ath-MIR5639-201", "noid"],
                "RNAcentral": ["URS000078C7F0"],
            },
            species="Arabidopsis thaliana",
            common_name="thale-cress",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "eudicotyledons; Gunneridae; Pentapetalae; rosids; malvids; "
                "Brassicales; Brassicaceae; Camelineae; Arabidopsis; "
                "Arabidopsis thaliana"
            ),
            gene="AT2G08015",
            locus_tag="ath-MIR5639",
            description=(
                "Arabidopsis thaliana (thale-cress) arabidopsis thaliana"
                " miR5639 stem-loop"
            ),
            references=[pubs.reference(29092050)],
        )
    )


# OBART09G00240.1  transposable_elements
# ENSRNA049475598-T1  sense_intronic
@pytest.mark.skip()
def test_skips_transposable_elements(oryza_9):
    pass


def test_can_parse_rice_trna(oryza_9):
    val = attr.asdict(helpers.entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049456349-T1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="ENSRNA049456349-T1",
            accession="ENSEMBL_PLANTS:ENSRNA049456349-T1",
            ncbi_tax_id=65489,
            database="ENSEMBL_PLANTS",
            sequence="TCCGTTGTAGTCTAGCTGGTTAGGATACTCGGCTCTCACCCGAGAGACCCGGGTTCGAGTCCCGGCAACGGAA",
            regions=[
                dat.SequenceRegion(
                    chromosome="9",
                    strand=1,
                    exons=[dat.Exon(start=31, stop=103)],
                    assembly_id="O.barthii_v1",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "GO": ["GO:0006412", "GO:0030533"],
                "RNAcentral": ["URS00001F06B0"],
                "TRNASCAN_SE_trans_name": ["tRNA-Glu-201"],
            },
            species="Oryza barthii",
            common_name="African wild rice",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; "
                "Oryzeae; Oryzinae; Oryza; Oryza barthii"
            ),
            gene="ENSRNA049456349",
            locus_tag="tRNA-Glu",
            description="Oryza barthii (African wild rice) misc RNA tRNA-Glu",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_rice_snorna(oryza_9):
    val = attr.asdict(helpers.entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049475670-T1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="ENSRNA049475670-T1",
            accession="ENSEMBL_PLANTS:ENSRNA049475670-T1",
            ncbi_tax_id=65489,
            database="ENSEMBL_PLANTS",
            sequence="AAAAAAGCAGGATGCTGTGTTCTCTATAAGCAGTGTCCTCGTAAATTTTAGGAACATGTTTCATCGTTATTGGGTGAACCGTTGGGCTATTCAATGTCCATTGGTTCAGTAAATGATGGCACATTT",
            regions=[
                dat.SequenceRegion(
                    chromosome="9",
                    strand=-1,
                    exons=[dat.Exon(start=164, stop=289)],
                    assembly_id="O.barthii_v1",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "GO": ["GO:0005730", "GO:0006396"],
                "RFAM_trans_name": ["snoR74-201"],
                "RNAcentral": ["URS0000653525"],
            },
            species="Oryza barthii",
            common_name="African wild rice",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; "
                "Oryzeae; Oryzinae; Oryza; Oryza barthii"
            ),
            gene="ENSRNA049475670",
            locus_tag="snoR74",
            description="Oryza barthii (African wild rice) misc RNA snoR74",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_rice_pre_mirna(oryza_9):
    val = attr.asdict(helpers.entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049475651-T1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="ENSRNA049475651-T1",
            accession="ENSEMBL_PLANTS:ENSRNA049475651-T1",
            ncbi_tax_id=65489,
            database="ENSEMBL_PLANTS",
            sequence="CCTCCCCGCCGGACCTCCCAGTGAGGAGGCTAGGGCCGCCAGGTCCGGTGATCCCATTCTCCTTGCCGGCGGATTCTGCGCCCTAGA",
            regions=[
                dat.SequenceRegion(
                    chromosome="9",
                    strand=-1,
                    exons=[dat.Exon(start=350, stop=436)],
                    assembly_id="O.barthii_v1",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "GO": ["GO:0035195"],
                "RFAM_trans_name": ["MIR1846-201"],
                "RNAcentral": ["URS0000BFFA82"],
            },
            species="Oryza barthii",
            common_name="African wild rice",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; "
                "Oryzeae; Oryzinae; Oryza; Oryza barthii"
            ),
            gene="ENSRNA049475651",
            locus_tag="MIR1846",
            description="Oryza barthii (African wild rice) misc RNA MIR1846",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_rice_u6(oryza_9):
    val = attr.asdict(helpers.entry_for(oryza_9, "ENSEMBL_PLANTS:ENSRNA049475710-T1"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="ENSRNA049475710-T1",
            accession="ENSEMBL_PLANTS:ENSRNA049475710-T1",
            ncbi_tax_id=65489,
            database="ENSEMBL_PLANTS",
            sequence="GTAGCTTATATACGCTGCTGTGCATAAAATTGAAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGAAGACGCACACAAATCGAGAAGTGGTCCAAATTTTT",
            regions=[
                dat.SequenceRegion(
                    chromosome="9",
                    strand=1,
                    exons=[dat.Exon(start=497, stop=603)],
                    assembly_id="O.barthii_v1",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={
                "GO": [
                    "GO:0000244",
                    "GO:0000353",
                    "GO:0005688",
                    "GO:0030621",
                    "GO:0046540",
                ],
                "RFAM_trans_name": ["U6-201"],
                "RNAcentral": ["URS0000C52597"],
            },
            species="Oryza barthii",
            common_name="African wild rice",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "Liliopsida; Poales; Poaceae; BOP clade; Oryzoideae; "
                "Oryzeae; Oryzinae; Oryza; Oryza barthii"
            ),
            gene="ENSRNA049475710",
            locus_tag="U6",
            description="Oryza barthii (African wild rice) misc RNA U6",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_barley_antisense(hordeum_pt):
    val = attr.asdict(
        helpers.entry_for(hordeum_pt, "ENSEMBL_PLANTS:ENSRNA049483195-T1")
    )
    assert val == attr.asdict(
        dat.Entry(
            primary_id="ENSRNA049483195-T1",
            accession="ENSEMBL_PLANTS:ENSRNA049483195-T1",
            ncbi_tax_id=112509,
            database="ENSEMBL_PLANTS",
            sequence="AATAACCAAATATAACACTGGGACTAAGGGTCAAATTGGTAATTTTTCTTACATCTCCCCCCCCAGGGGCCCAGGTATCATATACACCGCCAAAATAAAGAGCCTTGAGTACTAGAAGAAAAGCACCTAGACCTAACAAAATTAAGTGAATACCCAAAATTGTAGTCATTTTATTTCTATCTTTCCA",
            regions=[
                dat.SequenceRegion(
                    chromosome="Pt",
                    strand=-1,
                    exons=[dat.Exon(start=10076, stop=10262)],
                    assembly_id="IBSC_v2",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="antisense_RNA",
            url="",
            seq_version="1",
            species="Hordeum vulgare subsp. vulgare",
            common_name="two-rowed barley",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; Embryophyta; "
                "Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; "
                "Poales; Poaceae; BOP clade; Pooideae; Triticodae; "
                "Triticeae; Hordeinae; Hordeum; Hordeum vulgare subsp. vulgare"
            ),
            gene="ENSRNA049483195",
            locus_tag="IsrR",
            description="Hordeum vulgare subsp. vulgare (two-rowed barley) antisense RNA which regulates isiA expression",
            references=[pubs.reference(29092050)],
        )
    )


def test_can_parse_zea_lincrna(zea_7):
    # Zea mays moved from assembly B73_RefGen_v4 to Zm-B73-REFERENCE-NAM-5.0
    # (gene ids Zm00001d... -> Zm00001eb...) since this fixture was recorded;
    # re-pointed at a current misc_RNA gene from the new assembly.
    val = attr.asdict(helpers.entry_for(zea_7, "ENSEMBL_PLANTS:Zm00001eb302990_T001"))
    assert val == attr.asdict(
        dat.Entry(
            primary_id="Zm00001eb302990_T001",
            accession="ENSEMBL_PLANTS:Zm00001eb302990_T001",
            ncbi_tax_id=4577,
            database="ENSEMBL_PLANTS",
            sequence=(
                "ACAAATTTCAATTGTGCCAATTTTTGTGGGTATGCAAAGGAAAGGAGAAATAGAGGAAAT"
                "GATTATGGAGTCGAAGCATGCCATGGAGCTTTCATTGTATATTACATATGTTGGATGGA"
                "CGATAGAGGTGATAGGCAGTAAGCTAGTAGGTAGCTATGAATATTTGTTCTATAATTAA"
                "GTTTCTCACAGTTTTTGCATTTATGTTCTTGGATAAGAGTTGAAAAGATTGCTACTGCT"
                "ACTTCAGTTTTCAGTAGATTCTGCATAGAAGATCATTATGTTAATGTAGTGCAAATCC"
                "TTTGGGATTTAGGTGATAGAATTATATTTGCTTTGAAGATGCTAAAAATCTACCTGTA"
                "TTTCTTTCCTCCC"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="7",
                    strand=-1,
                    exons=[dat.Exon(start=31, stop=396)],
                    assembly_id="Zm-B73-REFERENCE-NAM-5.0",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                )
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            xref_data={"RNAcentral": ["URS00021F5853"]},
            species="Zea mays",
            common_name="maize",
            lineage=(
                "Eukaryota; Viridiplantae; Streptophyta; "
                "Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; "
                "Liliopsida; Poales; Poaceae; PACMAD clade; Panicoideae; "
                "Andropogonodae; Andropogoneae; Tripsacinae; Zea; Zea mays"
            ),
            gene="Zm00001eb302990",
            description="Zea mays (maize) misc RNA Zm00001eb302990",
            references=[pubs.reference(29092050)],
        )
    )


def test_does_not_generate_tair_for_others(zea_7):
    assert helpers.has_entry_for(zea_7, "TAIR:Zm00001eb302990_T001") is False
    assert helpers.has_entry_for(zea_7, "ENSEMBL_PLANTS:Zm00001eb302990_T001") is True


def test_does_not_create_ncRNA_rna_type_zea_7(zea_7):
    for entry in zea_7:
        assert entry.rna_type != "ncRNA"


def test_does_not_create_ncRNA_rna_type_cress_2(cress_2):
    for entry in cress_2:
        assert entry.rna_type != "ncRNA"


def test_does_not_create_ncRNA_rna_type_oryza_9(oryza_9):
    for entry in oryza_9:
        assert entry.rna_type != "ncRNA"


def test_does_not_create_ncRNA_rna_type_hordeum_pt(hordeum_pt):
    for entry in hordeum_pt:
        assert entry.rna_type != "ncRNA"
