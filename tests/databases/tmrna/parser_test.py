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


def test_can_parse_file(data):
    assert len(data) == 100


@pytest.mark.parametrize(
    "id,expected",
    [
        (
            "Paulinella__chromatophora.1:1",
            Entry(
                primary_id="Paulinella__chromatophora.1:1",
                accession="Paulinella__chromatophora.1:1",
                seq_version="1",
                ncbi_tax_id=39717,
                database="TMRNA_WEB",
                regions=[],
                sequence="GTTCGGTTATTGCCGAACTAGGTGGCTCACACCAATGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCAttattaGGGGCTGCAATGGTTTCGACGGGGCATCAGGAGGGTTACTGAAGCCTGCTCGGTAAGAGCAAATTAGTAACAgcgaacaacatcgttcgtttctcccgtcaagcggcccctgtggctgccTGACCCTAGATAGGGAGATGAGGTAAAGTCAGCCTTATAACCCAAATGACTCAAGGGGCCTGTAAGGGCCCCATCATTA",
                url="",
                rna_type="SO:0000584",
                features=[
                    SequenceFeature(
                        name="",
                        feature_type="",
                        location=[0, 70],
                        sequence="",
                        provider="TMRNA_WEBSITE",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="",
                        feature_type="",
                        location=[],
                        sequence="",
                        provider="TMRNA_WEBSITE",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="",
                        feature_type="",
                        location=[],
                        sequence="",
                        provider="TMRNA_WEBSITE",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="",
                        feature_type="",
                        location=[],
                        sequence="",
                        provider="TMRNA_WEBSITE",
                        metadata={},
                    ),
                    SequenceFeature(
                        name="",
                        feature_type="",
                        location=[],
                        sequence="",
                        provider="TMRNA_WEBSITE",
                        metadata={},
                    ),
                ],
            ),
        ),
        (
            "Paulinella__chromatophora.1:2",
            Entry(
                primary_id="Paulinella__chromatophora.1:2",
                accession="Paulinella__chromatophora.1:2",
                seq_version="1",
                ncbi_tax_id=39717,
                database="TMRNA_WEB",
                regions=[],
                sequence="GTTCGGTTATTGCCGAACTAGGTGGCTCACACCAATGTTTCGGACAGCGGTTCGATTCCGCTCAGCTCCAttattaGGGGCTGCAATGGTTTCGACGGGGCATCAGGAGGGTTACTGAAGCCTGCTCGGTAAGAGCAAATTAGTAACAgcgaacaacatcgttcgtttctcccgtcaagcggcccctgtggctgccTGACCCTAGATAGGGAGATGAGGTAAAGTCAGCCTTATAACCCAAATGACTCAAGGGGCCTGTAAGGGCCCCATCATTA",
                url="",
                rna_type="SO:0000584",
            ),
        ),
    ],
)
def test_can_parse_correctly(data, id, expected):
    assert data[id] == expected
