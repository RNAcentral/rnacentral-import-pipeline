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

from rnacentral_pipeline.databases.gencode import parser as gencode

from ..ensembl.helpers import parse_with_family as ensembl_parse
from . helpers import parse, entries_for, entry_for


@pytest.fixture(scope='module')  # pylint: disable=no-member
def human_12():
    return parse('data/ensembl/Homo_sapiens.GRCh38.90.chromosome.12.dat')


def known():
    with open('data/gencode/human.gff', 'rb') as raw:
        return gencode.gencode_transcripts(raw)


def test_can_load_all_known_transcripts(known):
    assert len(known) == 1000


def test_can_load_all_entries(known, human_12, family_file):
    assert len(list(gencode.parse(known, human_12, family_file))) == 100


@pytest.mark.parameterize('transcript,status', [  # pylint: disable=no-member
])
def test_knows_if_from_gencode(known, human_12, transcript_id, status):
    entry = ensembl_entry(human_12, transcript_id)
    assert gencode.from_gencode(known, entry) == status


def test_can_properly_transform_an_entry(known, human_12):
    assert entry_for(known, human_12, '') == dat.Entry()


# @pytest.mark.parameterize('transcript,status', [  # pylint: disable=no-member
#     ("ENST00000535849.1", True)
# ])
# def test_knows_if_entry_is_gencode(human_12, transcript_id, status):
#     assert bool(entries_for(human_12, transcript_id)) == status


# def test_it_uses_correct_primary_id(human_12):
#     entry = entry_for(human_12, "ENST00000400706.3")
#     assert entry.primary_id == 'ENST00000400706'


# def test_it_adds_xref_to_ensembl(human_12):
#     entry = entry_for(human_12, "ENST00000540226.1")
#     assert 'ENST00000540226.1' in entry.xref_data['Ensembl']


# def test_it_gets_all_gencode_entries(human_12):
#     assert len(human_12) == 2898


# @pytest.mark.parameterize('transcript,acc', [  # pylint: disable=no-member
#     ("ENST00000540868.1", 'GENCODE:ENST00000540868.1'),
# ])
# def test_it_sets_accession_correctly(human_12, transcript, acc):
#     assert entry_for(human_12, transcript).accession == acc
