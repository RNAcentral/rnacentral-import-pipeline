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

import pytest

from rnacentral_pipeline.databases.data import Entry, SequenceFeature
from rnacentral_pipeline.databases.tmrna import parser


@pytest.fixture(scope="module")
def data():
    with open("data/tmrna/example.tsv", "r") as raw:
        data = {}
        for item in parser.parse(raw):
            assert item.primary_id not in data
            data[item.primary_id] = item
        return data


@pytest.mark.network
def test_can_parse_file(data):
    assert len(data) == 108


@pytest.mark.network
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
