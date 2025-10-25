# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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

from unittest.mock import patch

import pytest

from rnacentral_pipeline.databases.data import Entry, SequenceFeature
from rnacentral_pipeline.databases.tmrna import parser

# Taxonomy mappings extracted from test data to avoid network calls
# Maps GTDB lineage strings to NCBI taxonomy IDs
TAXONOMY_MAPPINGS = {
    "Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae,Paulinella chromatophora": 39717,
    "Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae,Paulinella longichromatophora": 1708747,
    "Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae,Paulinella micropora": 1928728,
    "Plastids,d__Eukaryota,Schizocladia ischiensis": 196139,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae,Cryptomonas curvata": 233186,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae,Cryptomonas paramecium": 2898,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Cryptomonadales,f__Cryptomonadaceae,Cryptomonas pyrenoidifera": 233184,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Pyrenomonadales,f__Chroomonadaceae,Chroomonas placoidea": 173977,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Pyrenomonadales,f__Geminigeraceae,Guillardia theta": 55529,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Pyrenomonadales,f__Geminigeraceae,Teleaulax amphioxeia": 77931,
    "Plastids,d__Eukaryota,c__Cryptophyceae,o__Pyrenomonadales,f__Pyrenomonadaceae,Rhodomonas salina": 3034,
    "Plastids,d__Eukaryota,c__Glaucocystophyceae,f__Cyanophoraceae,Cyanophora biloba": 1489483,
    "Plastids,d__Eukaryota,c__Glaucocystophyceae,f__Cyanophoraceae,Cyanophora paradoxa": 2762,
    "Plastids,d__Eukaryota,c__Glaucocystophyceae,f__Cyanophoraceae,Cyanophora sudae": 1522369,
    "Plastids,d__Eukaryota,k__Alveolata,c__Dinophyceae,o__Peridiniales,f__Kryptoperidiniaceae,Durinskia baltica": 59809,
    "Plastids,d__Eukaryota,k__Alveolata,c__Dinophyceae,o__Peridiniales,f__Kryptoperidiniaceae,Kryptoperidinium foliaceum": 160619,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Bolidophyceae,o__Parmales,f__Triparmaceae,Bolidomonas sp": 722751,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Bolidophyceae,o__Parmales,f__Triparmaceae,Triparma laevis": 1534972,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Dictyochophyceae,o__Dictyochales,Dictyocha speculum": 3111310,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Dictyochophyceae,o__Florenciellales,Florenciella parvula": 236787,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Dictyochophyceae,o__Pedinellales,Pseudopedinella elastica": 35684,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Dictyochophyceae,o__Rhizochromulinales,Rhizochromulina marina": 1034831,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,Eustigmatophyceae sp": 5747,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Characiopsis acuta": 2040456,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Chlorobotrys sp": 2974601,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Lietzensia polymorpha": 2962110,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Neustupella aerophytica": 2962111,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Vischeria punctata": 643629,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Chlorobotryaceae,Vischeria sp": 2974601,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Monodopsis sp": 425072,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis gaditana": 72520,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis granulata": 43926,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis limnetica": 120807,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis oceanica": 145522,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis oculata": 43925,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Monodopsidaceae,Nannochloropsis salina": 2511165,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Eustigmatales,f__Neomonodaceae,Pseudellipsoidion edaphicum": 1431838,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Eustigmatophyceae,o__Goniochloridales,f__Goniochloridaceae,Trachydiscus minutus": 1032745,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Olisthodiscophyceae,f__Olisthodiscaceae,Olisthodiscus luteus": 83000,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Pelagophyceae,o__Pelagomonadales,Aureococcus anophagefferens": 44056,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Pelagophyceae,o__Pelagomonadales,Aureoumbra lagunensis": 44058,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Pelagophyceae,o__Pelagomonadales,Pelagomonas sp": 54409,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Dictyotales,f__Dictyotaceae,Dictyopteris divaricata": 156996,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Acinetosporaceae,Pylaiella littoralis": 2885,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Chordariaceae,Cladosiphon okamuranus": 309737,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Chordariaceae,Pleurocladia lacustris": 246121,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Ectocarpaceae,Ectocarpus siliculosus": 2880,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Ishigeaceae,Ishige okamurae": 233772,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Scytosiphonaceae,Colpomenia_sinuosa": 2891,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Scytosiphonaceae,Endarachne binghamiae": 698476,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Scytosiphonaceae,Scytosiphon canaliculatus": 2567908,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Scytosiphonaceae,Scytosiphon lomentaria": 27967,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Ectocarpales,f__Scytosiphonaceae,Scytosiphon promiscuus": 1403536,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Fucaceae,Fucus vesiculosus": 49266,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Fucaceae,Silvetia siliquosa": 93837,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Coccophora langsdorfii": 74099,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum confusum": 74091,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum fulvellum": 3016,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum fusiforme": 590727,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum horneri": 74089,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum macrocarpum": 74092,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum plagiophyllum": 1436148,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum polycystum": 127578,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum siliquastrum": 127572,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Fucales,f__Sargassaceae,Sargassum thunbergii": 127542,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Agaraceae,Costaria costata": 2872,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Alariaceae,Alaria crassifolia": 98220,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Alariaceae,Alaria crispa": 441892,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Alariaceae,Alaria marginata": 98221,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Alariaceae,Alaria praelonga": 88159,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Alariaceae,Undaria pinnatifida": 74381,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Laminaria digitata": 80365,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Laminaria rodriguezii": 1740620,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Laminaria solidungula": 309363,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Macrocystis integrifolia": 169774,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Macrocystis pyrifera": 35122,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Saccharina japonica": 88149,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Laminariaceae,Saccharina latissima": 309358,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Lessoniaceae,Ecklonia radiata": 309355,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Lessoniaceae,Lessonia flavicans": 169771,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Phaeophyceae,o__Laminariales,f__Lessoniaceae,Lessonia spicata": 1899210,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,Haramonas pauciplastida": 478668,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,o__Chattonellales,f__Chattonellaceae,Chattonella marina": 90936,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,o__Chattonellales,f__Chattonellaceae,Heterosigma akashiwo": 2829,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,o__Chattonellales,f__Vacuolariaceae,Gonyostomum semen": 375454,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,o__Chattonellales,f__Vacuolariaceae,Merotricha bacillata": 658122,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Raphidophyceae,o__Chattonellales,f__Vacuolariaceae,Vacuolaria virescens": 44451,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Synurophyceae,o__Synurales,f__Mallomonadaceae,Mallomonas splendens": 52552,
    "Plastids,d__Eukaryota,k__Ochrophyta,c__Synurophyceae,o__Synurales,f__Mallomonadaceae,Synura uvella": 52557,
    "Plastids,d__Eukaryota,k__Ochrophyta,p__Bacillariophyta,c__Bacillariophyceae,f__Entomoneidaceae,Entomoneis sp": 186041,
    "d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095,WVXT01 sp009619095": 1869227,
    "d__Bacteria,p__Verrucomicrobiota,c__Verrucomicrobiae,o__Opitutales,f__UBA953,UBA953__sp003569205,UBA953 sp003569205": 415000,
}


@pytest.fixture(scope="module")
def mock_taxonomy():
    """Mock the taxonomy lookup to avoid network calls to ENA/UniProt APIs"""
    def mock_phylogeny_to_taxid(lineage: str) -> int:
        if lineage in TAXONOMY_MAPPINGS:
            return TAXONOMY_MAPPINGS[lineage]
        # Raise error for unmapped lineages - should not happen with complete test data
        raise ValueError(f"Unmapped lineage in test data: {lineage}")

    with patch("rnacentral_pipeline.databases.helpers.gtdb.phylogeny_to_taxid", side_effect=mock_phylogeny_to_taxid):
        yield


@pytest.fixture(scope="module")
def data(mock_taxonomy):
    with open("data/tmrna/example.tsv", "r") as raw:
        data = {}
        for item in parser.parse(raw):
            assert item.primary_id not in data
            data[item.primary_id] = item
        return data


def test_can_parse_file(data):
    assert len(data) == 108


@pytest.mark.parametrize(
    "id,expected",
    [
        (
            "tmrna:CP000815.1/744167-744441",
            Entry(
                primary_id="tmrna:CP000815.1/744167-744441",
                accession="tmrna:CP000815.1/744167-744441",
                seq_version="1",
                ncbi_tax_id=39717,
                database="TMRNA_WEB",
                regions=[],
                sequence="GTTCGGTTATTGCCGAACTAGGTGGCTCACACCAATGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCATTATTAGGGGCTGCAATGGTTTCGACGGGGCATCAGGAGGGTTACTGAAGCCTGCTCGGTAAGAGCAAATTAGTAACAGCGAACAACATCGTTCGTTTCTCCCGTCAAGCGGCCCCTGTGGCTGCCTGACCCTAGATAGGGAGATGAGGTAAAGTCAGCCTTATAACCCAAATGACTCAAGGGGCCTGTAAGGGCCCCATCATTA",
                url="",
                rna_type="SO:0000584",
                inference="Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae",
                parent_accession="CP000815.1",
                note_data={
                    "tmrna_form": "Permuted",
                },
                description="Paulinella chromatophora Permuted tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_acceptor",
                        feature_type="tmrna_acceptor",
                        location=[0, 70],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_ivs",
                        feature_type="tmrna_ivs",
                        location=[70, 76],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_coding_region",
                        feature_type="tmrna_coding_region",
                        location=[76, 275],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[148, 199],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "Coding sequence",
                            "has_stop": True,
                            "coding_sequence": "ANNIVRFSRQAAPVAA",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[70, 73],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                ],
            ),
        ),
        (
            "tmrna:WJOT01000061.1/9381-9976",
            Entry(
                primary_id="tmrna:WJOT01000061.1/9381-9976",
                accession="tmrna:WJOT01000061.1/9381-9976",
                seq_version="1",
                ncbi_tax_id=1869227,
                database="TMRNA_WEB",
                regions=[],
                sequence="GGGGGCGATCAGGTTTCGACAGGAATAAAGGAGGCAAGGACGGCAGGTCGAGGTTTGTCGAAGGCTCGTTAATCAATCGACAACAAAAACTAAGTGCTGACACTAAGTTAGCATTAGCCGCGTAAGCGGACACGCTCACCTCTTTTTGCCCATTGGATGGGGATGAGTGTCAGAGAGATGGGATAGTTCCGGCTTTTGCCTTGGAAGTCGGGATGAGATTTAAAGGCTGCTCCGTCTCGAATCCTGTCTTTGGGTAATAAGGCGGGGGAGATTCCAAACAAAGACTAAACCTGTAGATGTCCTGCTGAAATATTTCTGGACGCGGGTTCGATCAAAAGGTCGACTTCCGAAGTAATTCGGGATGCAAAACCTGGCTTACATCAGGGAAACCCTAATCCTGATAATACAGGAAAGGGTAATCCTGAGGGGCGAAAGCCCTGCAGAGACTGTAATGATAGGAGGTGATATCGTTGGACGAAGTTGAACCAACGATAAATCTATGCCTCCTATAACACGTCAGGCTCCCCGATTGGGGATGAAGAGATAGTCCATACCACTTAGAAATAAATGGAGTATATGTCCCGCCGCCTCCACCA",
                url="",
                rna_type="SO:0000584",
                inference="d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095",
                parent_accession="WJOT01000061.1",
                note_data={
                    "tmrna_form": "GpIintron",
                },
                description="WVXT01 sp009619095 GpIintron tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_exon",
                        feature_type="tmrna_exon",
                        location=[0, 332],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[95, 125],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "Coding sequence",
                            "has_stop": True,
                            "coding_sequence": "ADTKLALAA",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_gpi",
                        feature_type="tmrna_gpi",
                        location=[332, 579],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_exon",
                        feature_type="tmrna_exon",
                        location=[579, 593],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[593, 596],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                ],
            ),
        ),
        (
            "tmrna:BGOK01000023.1/25516-25150",
            Entry(
                primary_id="tmrna:BGOK01000023.1/25516-25150",
                accession="tmrna:BGOK01000023.1/25516-25150",
                seq_version="1",
                ncbi_tax_id=415000,
                database="TMRNA_WEB",
                regions=[],
                sequence="GGGGGTGTTTCCGGATTCGATTCCCATCAGATATCATGACGGCATGCAGAGGACGTCAGCTCCCTCTTAAATCCAGCTGGCAGCATATAACTGCTAAAAACAACCATCACGTTCGAAAGCGACAACGCTTTCCGCACTCGCTGCTTAATTACAGCGACGGGCCCAGGACGCGCCTGCATTCCTGGATCAGCTCGGTTAACGTGGGCGATCCCCTCCTCTAGTTTCTATCTCCGGGGTAAGTGAGAATGAAGATAGGCCGCACGGGTGCGCATTCCCAAACGCGGTCGAGATCAATAATGACGCTAAGCATGTAGAAGATGTGACGTAAGGATAGGAAGACGCGGGTCGACTCCCGCCACCTCCACCA",
                url="",
                rna_type="SO:0000584",
                inference="d__Bacteria,p__Verrucomicrobiota,c__Verrucomicrobiae,o__Opitutales,f__UBA953,UBA953__sp003569205",
                parent_accession="BGOK01000023.1",
                note_data={
                    "tmrna_form": "Standard",
                },
                description="UBA953 sp003569205 Standard tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_body",
                        feature_type="tmrna_body",
                        location=[0, 364],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[92, 148],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "frameshift",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[364, 367],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                ],
            ),
        ),
    ],
)
def test_can_parse_correctly(data, id, expected):
    assert data[id] == expected
