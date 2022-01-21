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

from rnacentral_pipeline.databases.ensembl import helpers

from .helpers import features, feature_for, first_feature_for


@pytest.fixture(scope="module")  # pylint: disable=no-member
def human_12():
    return features("data/ensembl/Homo_sapiens.GRCh38.chromosome.12.dat")


def test_it_gets_transcript_id(human_12):
    assert (
        helpers.transcript(feature_for(human_12, "ENST00000535849.1"))
        == "ENST00000535849.1"
    )


@pytest.mark.parametrize(
    "transcript_id,notes",
    [
        # ('ENST00000540907.11', [
        #     "processed_transcript",
        #     "transcript_id=ENST00000544511.1"
        # ]),
        ("ENST00000540226.1", ["antisense_RNA"])
    ],
)
def test_it_can_get_notes(human_12, transcript_id, notes):
    assert helpers.notes(feature_for(human_12, transcript_id)) == notes


@pytest.mark.parametrize(
    "transcript_id,note_data",
    [
        # ('ENST00000540907.11', {'transcript_id': ['ENST00000540907.11']}),
        ("ENST00000544511.1", {"transcript_id": ["ENST00000544511.1"]})
    ],
)
def test_it_can_get_grouped_notes(human_12, transcript_id, note_data):
    feature = feature_for(human_12, transcript_id)
    assert helpers.note_data(feature) == note_data


@pytest.mark.parametrize(
    "transcript_id,status",
    [
        ("ENST00000535572.5", False),
        ("ENST00000408512.1", True),
        ("ENST00000358495.7", False),
        ("ENST00000535376.5", False),
        ("ENST00000623153.1", False),
        ("ENST00000546223.1", True),
        ("ENST00000358495.7", False),
        ("ENST00000358495.7", False),
        ("ENST00000481052.5", False),
        ("ENST00000430095.6", False),
        ("ENST00000364606.1", True),
    ],
)
def test_can_detect_if_is_noncoding(human_12, transcript_id, status):
    feature = first_feature_for(human_12, transcript_id)
    assert helpers.is_ncrna(feature) == status
