# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

import os

import pytest

from rnacentral_pipeline.databases import data as dt
from rnacentral_pipeline.databases.ribovision import parser


def test_produces_expected_mapping():
    with open("data/ribovision/example.html", "r") as raw:
        data = parser.extract_mapping(raw)
    assert data == {
        "4V9D": (
            "E. coli",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        ),
        "1VY4": (
            "T. thermophilus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#THET8_LSU&THET8_SSU",
        ),
        "5NJT": (
            "B. subtilis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#BACSU_LSU",
        ),
        "4IOA": (
            "D. radiodurance",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DEIRA_LSU",
        ),
        "5O60": (
            "M. smegmatis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#MYCSM_LSU",
        ),
        "5V7Q": (
            "M. tuberculosis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#MYCTU_LSU",
        ),
        "4WF9": (
            "S. aureus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#STAAU_LSU",
        ),
        "4V6F": (
            "H. marismortui",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HALMA_LSU",
        ),
        "4V6U": (
            "P. furiousus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#PYRFU_LSU&PYRFU_SSU",
        ),
        "4V88": (
            "S. cerevisiae",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#YEAST_LSU&YEAST_SSU",
        ),
        "4V8P": (
            "T. thermophylia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#TETTH_LSU",
        ),
        "4V6W": (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        "4V6X": (
            "H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_LSU&HUMAN_SSU",
        ),
        "3J7Y": (
            "mt-H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_MTLSU",
        ),
        "6ERI": (
            "cl-S. oleracia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#SPIOL_LSU&SPIOL_SSU",
        ),
    }


@pytest.mark.db
def produces_expected_data():
    with open("data/ribovision/example.html", "r") as raw:
        data = next(parser.extract_mapping(raw))
    assert data == dt.Entry(
        primary_id="ribovision:4V9D",
        accession="ribovision:4V9D",
        ncbi_tax_id=562,
        database="RIBOVISION",
        sequence="",
        regions=[],
        rna_type="",
        url="http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        seq_version="1",
        description="E. coli rRNA",
    )
