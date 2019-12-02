# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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
import tempfile
from contextlib import contextmanager

import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.genecards_suite.genecards import parse
from rnacentral_pipeline.databases.genecards_suite.core import lookup 


@contextmanager
def known(handle):
    with tempfile.NamedTemporaryFile() as tmp:
        lookup.write(handle, os.environ['PGDATABASE'], 'URSid', tmp)
        tmp.seek(0)
        yield tmp


@pytest.fixture(scope='module')
def simple_data():
    with open('data/genecards/data.tsv', 'r') as raw:
        with known(raw) as indexed:
            raw.seek(0)
            return list(parse(raw, indexed))


def test_can_parse_all_entries(simple_data):
    assert len(simple_data) == 100


def test_can_create_unique_primary_ids(simple_data):
    data = [d.primary_id for d in simple_data]
    assert len(data) == 100


def test_can_create_unique_accessions(simple_data):
    data = [d.accession for d in simple_data]
    assert len(data) == 100


def test_can_create_correct_data(simple_data):
    assert simple_data[0] == data.Entry(
        primary_id='GENECARDS:1A9N_Q-015:URS00001EE9F1_9606',
        accession='GENECARDS:1A9N_Q-015:URS00001EE9F1_9606',
        ncbi_tax_id=9606,
        database='GENECARDS',
        sequence='ATTGCAGTACCTCCAGGAACGGTGCAC',
        regions=[],
        rna_type='misc_RNA',
        url='https://www.genecards.org/cgi-bin/carddisp.pl?gene=1A9N_Q-015',
        seq_version=1,
        gene='1A9N_Q-015',
        description='Homo sapiens (human) miscellaneous RNA',
    )
