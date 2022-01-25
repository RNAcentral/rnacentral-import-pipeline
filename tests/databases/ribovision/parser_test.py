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
        ("4V9D", "CA"): (
            "E. coli",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        ),
        ("4V9D", "CB"): (
            "E. coli",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        ),
        ("4V9D", "BA"): (
            "E. coli",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        ),
        ("1VY4", "Z"): (
            "T. thermophilus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#THET8_LSU&THET8_SSU",
        ),
        ("1VY4", "EB"): (
            "T. thermophilus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#THET8_LSU&THET8_SSU",
        ),
        ("1VY4", "AA"): (
            "T. thermophilus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#THET8_LSU&THET8_SSU",
        ),
        ("5NJT", "U"): (
            "B. subtilis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#BACSU_LSU",
        ),
        ("4IOA", "X"): (
            "D. radiodurance",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DEIRA_LSU",
        ),
        ("5O60", "A"): (
            "M. smegmatis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#MYCSM_LSU",
        ),
        ("5V7Q", "A"): (
            "M. tuberculosis",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#MYCTU_LSU",
        ),
        ("4WF9", "X"): (
            "S. aureus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#STAAU_LSU",
        ),
        ("1JJ2", "0"): (
            "H. marismortui",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HALMA_LSU",
        ),
        ("1JJ2", "9"): (
            "H. marismortui",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HALMA_LSU",
        ),
        ("4V6U", "A2"): (
            "P. furiousus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#PYRFU_LSU&PYRFU_SSU",
        ),
        ("4V6U", "B1"): (
            "P. furiousus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#PYRFU_LSU&PYRFU_SSU",
        ),
        ("4V6U", "B3"): (
            "P. furiousus",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#PYRFU_LSU&PYRFU_SSU",
        ),
        ("4V88", "A4"): (
            "S. cerevisiae",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#YEAST_LSU&YEAST_SSU",
        ),
        ("4V88", "A1"): (
            "S. cerevisiae",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#YEAST_LSU&YEAST_SSU",
        ),
        ("4V88", "A3"): (
            "S. cerevisiae",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#YEAST_LSU&YEAST_SSU",
        ),
        ("4V88", "A2"): (
            "S. cerevisiae",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#YEAST_LSU&YEAST_SSU",
        ),
        ("4V8P", "B2"): (
            "T. thermophylia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#TETTH_LSU",
        ),
        ("4V8P", "A1"): (
            "T. thermophylia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#TETTH_LSU",
        ),
        ("4V8P", "B3"): (
            "T. thermophylia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#TETTH_LSU",
        ),
        ("4V6W", "A8"): (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        ("4V6W", "A9"): (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        ("4V6W", "A5"): (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        ("4V6W", "A7"): (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        ("4V6W", "B2"): (
            "D. melanogaster",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#DROME_LSU&DROME_SSU",
        ),
        ("4V6X", "A8"): (
            "H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_LSU&HUMAN_SSU",
        ),
        ("4V6X", "A5"): (
            "H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_LSU&HUMAN_SSU",
        ),
        ("4V6X", "A7"): (
            "H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_LSU&HUMAN_SSU",
        ),
        ("4V6X", "B2"): (
            "H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_LSU&HUMAN_SSU",
        ),
        ("3J7Y", "A"): (
            "mt-H. sapiens",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#HUMAN_MTLSU",
        ),
        ("6ERI", "AA"): (
            "cl-S. oleracia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#SPIOL_LSU&SPIOL_SSU",
        ),
        ("6ERI", "AX"): (
            "cl-S. oleracia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#SPIOL_LSU&SPIOL_SSU",
        ),
        ("6ERI", "BA"): (
            "cl-S. oleracia",
            "http://apollo.chemistry.gatech.edu/RiboVision2/#SPIOL_LSU&SPIOL_SSU",
        ),
    }


@pytest.mark.db
def test_produces_all_data():
    with open("data/ribovision/example.html", "r") as raw:
        data = list(parser.parse(raw, os.environ["PGDATABASE"]))
    assert len(data) == 36


@pytest.mark.db
def test_produces_expected_data():
    with open("data/ribovision/example.html", "r") as raw:
        data = list(parser.parse(raw, os.environ["PGDATABASE"]))
        data = data[0]
    assert data == dt.Entry(
        primary_id="ribovision:4V9D_CA",
        accession="ribovision:4V9D_CA",
        ncbi_tax_id=562,
        database="RIBOVISION",
        sequence="GGTTAAGCGACTAAGCGTACACGGTGGATGCCCTGGCAGTCAGAGGCGATGAAGGACGTGCTAATCTGCGATAAGCGTCGGTAAGGTGATATGAACCGTTATAACCGGCGATTTCCGAATGGGGAAACCCAGTGTGTTTCGACACACTATCATTAACTGAATCCATAGGTTAATGAGGCGAACCGGGGGAACTGAAACATCTAAGTACCCCGAGGAAAAGAAATCAACCGAGATTCCCCCAGTAGCGGCGAGCGAACGGGGAGCAGCCCAGAGCCTGAATCAGTGTGTGTGTTAGTGGAAGCGTCTGGAAAGGCGCGCGATACAGGGTGACAGCCCCGTACACAAAAATGCACATGCTGTGAGCTCGATGAGTAGGGCGGGACACGTGGTATCCTGTCTGAATATGGGGGGACCATCCTCCAAGGCTAAATACTCCTGACTGACCGATAGTGAACCAGTACCGTGAGGGAAAGGCGAAAAGAACCCCGGCGAGGGGAGTGAAAAAGAACCTGAAACCGTGTACGTACAAGCAGTGGGAGCACGCTTAGGCGTGTGACTGCGTACCTTTTGTATAATGGGTCAGCGACTTATATTCTGTAGCAAGGTTAACCGAATAGGGGAGCCGAAGGGAAACCGAGTCTTAACTGGGCGTTAAGTTGCAGGGTATAGACCCGAAACCCGGTGATCTAGCCATGGGCAGGTTGAAGGTTGGGTAACACTAACTGGAGGACCGAACCGACTAATGTTGAAAAATTAGCGGATGACTTGTGGCTGGGGGTGAAAGGCCAATCAAACCGGGAGATAGCTGGTTCTCCCCGAAAGCTATTTAGGTAGCGCCTCGTGAATTCATCTCCGGGGGTAGAGCACTGTTTCGGCAAGGGGGTCATCCCGACTTACCAACCCGATGCAAACTGCGAATACCGGAGAATGTTATCACGGGAGACACACGGCGGGTGCTAACGTCCGTCGTGAAGAGGGAAACAACCCAGACCGCCAGCTAAGGTCCCAAAGTCATGGTTAAGTGGGAAACGATGTGGGAAGGCCCAGACAGCCAGGATGTTGGCTTAGAAGCAGCCATCATTTAAAGAAAGCGTAATAGCTCACTGGTCGAGTCGGCCTGCGCGGAAGATGTAACGGGGCTAAACCATGCACCGAAGCTGCGGCAGCGACGCTTATGCGTTGTTGGGTAGGGGAGCGTTCTGTAAGCCTGCGAAGGTGTGCTGTGAGGCATGCTGGAGGTATCAGAAGTGCGAATGCTGACATAAGTAACGATAAAGCGGGTGAAAAGCCCGCTCGCCGGAAGACCAAGGGTTCCTGTCCAACGTTAATCGGGGCAGGGTGAGTCGACCCCTAAGGCGAGGCCGAAAGGCGTAGTCGATGGGAAACAGGTTAATATTCCTGTACTTGGTGTTACTGCGAAGGGGGGACGGAGAAGGCTATGTTGGCCGGGCGACGGTTGTCCCGGTTTAAGCGTGTAGGCTGGTTTTCCAGGCAAATCCGGAAAATCAAGGCTGAGGCGTGATGACGAGGCACTACGGTGCTGAAGCAACAAATGCCCTGCTTCCAGGAAAAGCCTCTAAGCATCAGGTAACATCAAATCGTACCCCAAACCGACACAGGTGGTCAGGTAGAGAATACCAAGGCGCTTGAGAGAACTCGGGTGAAGGAACTAGGCAAAATGGTGCCGTAACTTCGGGAGAAGGCACGCTGATATGTAGGTGAGGTCCCTCGCGGATGGAGCTGAAATCAGTCGAAGATACCAGCTGGCTGCAACTGTTTATTAAAAACACAGCACTGTGCAAACACGAAAGTGGACGTATACGGTGTGACGCCTGCCCGGTGCCGGAAGGTTAATTGATGGGGTTAGCGCAAGCGAAGCTCTTGATCGAAGCCCCGGTAAACGGCGGCCGTAACTATAACGGTCCTAAGGTAGCGAAATTCCTTGTCGGGTAAGTTCCGACCTGCACGAATGGCGTAATGATGGCCAGGCTGTCTCCACCCGAGACTCAGTGAAATTGAACTCGCTGTGAAGATGCAGTGTACCCGCGGCAAGACGGAAAGACCCCGTGAACCTTTACTATAGCTTGACACTGAACATTGAGCCTTGATGTGTAGGATAGGTGGGAGGCTTTGAAGTGTGGACGCCAGTCTGCATGGAGCCGACCTTGAAATACCACCCTTTAATGTTTGATGTTCTAACGTTGACCCGTAATCCGGGTTGCGGACAGTGTCTGGTGGGTAGTTTGACTGGGGCGGTCTCCTCCTAAAGAGTAACGGAGGAGCACGAAGGTTGGCTAATCCTGGTCGGACATCAGGAGGTTAGTGCAATGGCATAAGCCAGCTTGACTGCGAGCGTGACGGCGCGAGCAGGTGCGAAAGCAGGTCATAGTGATCCGGTGGTTCTGAATGGAAGGGCCATCGCTCAACGGATAAAAGGTACTCCGGGGATAACAGGCTGATACCGCCCAAGAGTTCATATCGACGGCGGTGTTTGGCACCTCGATGTCGGCTCATCACATCCTGGGGCTGAAGTAGGTCCCAAGGGTATGGCTGTTCGCCATTTAAAGTGGTACGCGAGCTGGGTTTAGAACGTCGTGAGACAGTTCGGTCCCTATCTGCCGTGGGCGCTGGAGAACTGAGGGGGGCTGCTCCTAGTACGAGAGGACCGGAGTGGACGCATCACTGGTGTTCGGGTTGTCATGCCAATGGCACTGCCCGGTAGCTAAATGCGGAAGAGATAAGTGCTGAAAGCATCTAAGCACGAAACTTGCCCCGAGATGAGTTCTCCCTGACCCTTTAAGGGTCCTGAAGGAACGTTGAAGACGACGACGTTGATAGGCCGGGTGTGTAAGCGCAGCGATGCGTTGAGCTAACCGGTACTAATGAACCGTGAGGCTTAACCT",
        regions=[],
        rna_type="SO:0001001",
        url="http://apollo.chemistry.gatech.edu/RiboVision2/#ECOLI_LSU&ECOLI_SSU",
        seq_version="1",
        description="E. coli rRNA",
    )
