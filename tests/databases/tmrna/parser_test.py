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

# Apply tmrna marker to all tests in this module
pytestmark = pytest.mark.tmrna

# TODO: Add validation test for Tag field format
# Issue found: Some source data has malformed Tag fields (e.g., "...AAALV11" instead of "...AAALV*")
# Should add test to validate:
#   1. Tag field ends with "*" (stop codon marker)
#   2. Tag field contains only valid amino acid codes + "*"
#   3. Parser handles malformed tags gracefully

# Taxonomy mappings extracted from test data to avoid network calls
# Maps Bacteriophage lineage strings to NCBI taxonomy IDs
# Updated for Bacteriophage-only dataset (25 entries, 8 unique lineages)
# All Bacteriophage entries resolve to tax_id 38018 ("Bacteriophages")
TAXONOMY_MAPPINGS = {
    "Bacteriophages,Phage Bacillus": 38018,
    "Bacteriophages,Phage Bordetella": 38018,
    "Bacteriophages,Phage Caulobacter": 38018,
    "Bacteriophages,Phage Cellulophaga": 38018,
    "Bacteriophages,Phage Enterococcus": 38018,
    "Bacteriophages,Phage Flavobacterium": 38018,
    "Bacteriophages,Phage Gordonia": 38018,
    "Bacteriophages,Phage Lactococcus": 38018,
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
    assert len(data) == 30  # 25 unique sequences, some with Count > 1 expand to multiple entries


@pytest.mark.parametrize(
    "id,expected",
    [
        (
            "tmrna:NC_023719.1/412060-412371",
            Entry(
                primary_id="tmrna:NC_023719.1/412060-412371",
                accession="tmrna:NC_023719.1/412060-412371",
                seq_version="1",
                ncbi_tax_id=38018,
                database="TMRNA_WEB",
                regions=[],
                sequence="GGGATTGTTAAGGATTCGACAGGGGAAAATAGTATAAGACATGATACTCGTAAAGCATACGTTACATGCTAAGTTAAATATAACTAACAACGAATTACAAGTAGCCTAGTTAATAGGTGCCAAAACAATAGAGTTGCTCTAACATCTATTTGAGGTTCAAAGATTTGAGCTATTTCTTACTTTAGAATAAAGTAAGTGGTGGAGCCTAGAGAAAACGCTGTAGAAACTAACTATTCAAATTGTAAATGGTATCAAACAAGTGTCTATATGAAAGTCTTCTGGACGCGGGTTCGACTCCCGCCAGTTCCATAT",
                url="",
                rna_type="SO:0000584",
                inference="Bacteriophages",
                parent_accession="NC_023719.1",
                note_data={
                    "tmrna_form": "Standard",
                },
                description="Phage Bacillus Standard tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_body",
                        feature_type="tmrna_body",
                        location=[0, 309],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[67, 109],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "Coding sequence",
                            "has_stop": True,
                            "coding_sequence": "AKLNITNNELQVA",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[309, 312],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                ],
            ),
        ),
        (
            "tmrna:NC_047861.1/46086-46387",
            Entry(
                primary_id="tmrna:NC_047861.1/46086-46387",
                accession="tmrna:NC_047861.1/46086-46387",
                seq_version="1",
                ncbi_tax_id=38018,
                database="TMRNA_WEB",
                regions=[],
                sequence="GGCCTCGGCCCGCGAGACGTGAAGTGAATGGCGTTCGGACGCGGGTTCAATTCCCGCCGGCTCCACCACTGGCAGCATTGGCGTCGAGCTCAGCATCCCCGGGGTCCGGATCACGGGTGAGTTCCAGTGCTGCCAGTCATGGGGCCGACCGGTTTCGACGGGCGCGCTGCAGCGGATCGGACTACTCGGCAAGCGGAGCCGTAAACCAAGCAAAATCGTAGACGCCAACGACGACTACATCCGCGCCGCTGCTTAAGCGGTGACGGCCTCCCCCGGCAGGTGGCAACAGAAGCCGGGGACCA",
                url="",
                rna_type="SO:0000584",
                inference="Bacteriophages",
                parent_accession="NC_047861.1",
                note_data={
                    "tmrna_form": "Permuted(beta)",
                },
                description="Phage Bordetella Permuted(beta) tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_acceptor",
                        feature_type="tmrna_acceptor",
                        location=[0, 65],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_ivs",
                        feature_type="tmrna_ivs",
                        location=[65, 140],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_coding_region",
                        feature_type="tmrna_coding_region",
                        location=[140, 302],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[223, 256],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "Coding sequence",
                            "has_stop": True,
                            "coding_sequence": "ANDDYIRAAA",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[65, 68],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                ],
            ),
        ),
        (
            "tmrna:NC_048047.1/149044-149586",
            Entry(
                primary_id="tmrna:NC_048047.1/149044-149586",
                accession="tmrna:NC_048047.1/149044-149586",
                seq_version="1",
                ncbi_tax_id=38018,
                database="TMRNA_WEB",
                regions=[],
                sequence="GGGACCTGGATGACGGAATATACGGAGACTGAACTGGACCGGGGTGCGATTCCCCGCGCCTCCACCACCAAGTCTCCCGGCCGGTCCCCACCCAGGGACCTGGGCGTCGCCTTCAGCGACGGCGCCCGATCTGAACGGCGCGGGAGACTTTGTTCTGGGGGCGATCAGCATCGACAGGCAGGGTAAGTAGAAGACGCGACACGGCATGGCTGCCGGTTCCTTCGGGAGCGCTCGGCCCACTTTGACAAATGTCAACGCTGCAAACGACAACGACGTCGTGGTTTCGACCATGACCTCGGTCAAGCTGGCCGCTTAACGGCTCGGGGTTCAGGAGGCGCCTTTTCACCCAAGGTCTCCAGCTTGCCCTAAGGGGCTTAGAGGACAGGCGGCGTAAGCCGCGCGATCAACACCAGGGGGCGGCCTGGGTCCCGAAACCGCCCACCCTCCTGCTTTCGCAGGTGGGATCAGAAGGTCGGCCTCAAAATAGGCGATCTGTCGGAAGAAACCCCGTTTCGCCGGCTTTCTGATCCCATCTACGAAAGC",
                url="",
                rna_type="SO:0000584",
                inference="Bacteriophages",
                parent_accession="NC_048047.1",
                note_data={
                    "tmrna_form": "Permuted(alpha)",
                },
                description="Phage Caulobacter Permuted(alpha) tmRNA",
                features=[
                    SequenceFeature(
                        name="tmrna_acceptor",
                        feature_type="tmrna_acceptor",
                        location=[0, 64],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_ivs",
                        feature_type="tmrna_ivs",
                        location=[64, 156],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_coding_region",
                        feature_type="tmrna_coding_region",
                        location=[156, 543],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="tmrna_tagcds",
                        feature_type="tmrna_tagcds",
                        location=[250, 316],
                        sequence="",
                        provider="TMRNA_WEB",
                        metadata={
                            "orf_summary": "Coding sequence",
                            "has_stop": True,
                            "coding_sequence": "VNAANDNDVVVSTMTSVKLAA",
                        },
                    ),
                    SequenceFeature(
                        name="tmrna_ccaequiv",
                        feature_type="tmrna_ccaequiv",
                        location=[64, 67],
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
