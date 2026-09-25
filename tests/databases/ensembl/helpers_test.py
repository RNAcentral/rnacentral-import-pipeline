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

import pytest

from rnacentral_pipeline.databases.ensembl.vertebrates import helpers

from .helpers import feature_for, features, first_feature_for


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_12():
    return features("test-data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat")


def test_it_gets_transcript_id(human_12):
    assert (
        helpers.transcript(feature_for(human_12, "ENST00000546223.1"))
        == "ENST00000546223.1"
    )


@pytest.mark.parametrize(
    "transcript_id,notes",
    [
        ("ENST00000540226.2", ["lncRNA"]),
    ],
)
def test_it_can_get_notes(human_12, transcript_id, notes):
    assert helpers.notes(feature_for(human_12, transcript_id)) == notes


@pytest.mark.parametrize(
    "transcript_id,note_data",
    [
        ("ENST00000432994.2", {"transcript_id": ["ENST00000432994.2"]}),
    ],
)
def test_it_can_get_grouped_notes(human_12, transcript_id, note_data):
    feature = feature_for(human_12, transcript_id)
    assert helpers.note_data(feature) == note_data


@pytest.mark.parametrize(
    "transcript_id,status",
    [
        ("ENST00000534526.7", False),  # mRNA, protein-coding
        ("ENST00000550091.5", True),  # misc_RNA, protein_coding_CDS_not_defined
        ("ENST00000494275.5", False),  # misc_RNA, retained_intron
        ("ENST00000531134.8", False),  # mRNA, protein-coding
        ("ENST00000611210.1", True),  # misc_RNA, misc_RNA
        ("ENST00000516089.1", True),  # misc_RNA, scaRNA
        ("ENST00000472289.5", False),  # mRNA, protein-coding
        ("ENST00000497153.5", True),  # misc_RNA, protein_coding_CDS_not_defined
    ],
)
def test_can_detect_if_is_noncoding(human_12, transcript_id, status):
    feature = first_feature_for(human_12, transcript_id)
    assert helpers.is_ncrna(feature) == status
