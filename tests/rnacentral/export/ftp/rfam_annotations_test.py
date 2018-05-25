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

import tempfile

import pytest
import psycopg2
import psycopg2.extras

from tasks.config import db

from rnacentral.export.ftp import rfam_annotations as rf

from tests.tasks.helpers import count


@pytest.mark.slowtest
def test_query_fetches_correct_data():
    connection = psycopg2.connect(db().pgloader_url())
    query = 'select * from (%s) t limit 1' % rf.QUERY
    cursor = connection.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute(query)
    assert dict(cursor.fetchone()) == {
        'upi': 'URS0000000001',
        'rfam_model_id': 'RF00177',
        'score': 109.4,
        'e_value': 2.1e-33,
        'sequence_start': 2,
        'sequence_stop': 200,
        'model_start': 29,
        'model_stop': 230,
        'long_name': 'Bacterial small subunit ribosomal RNA',
    }


@pytest.mark.slowtest
def test_query_produces_correct_counts():
    assert count(rf.QUERY) == 9683844  # 25403830


@pytest.mark.slowtest
def test_can_write_correct_file():
    connection = psycopg2.connect(db().pgloader_url())
    with tempfile.NamedTemporaryFile() as out:
        rf.write(connection, out)

        with open(out.name, 'r') as out:
            assert next(out) == 'URS0000000001\tRF00177\t109.4\t2.1e-33\t2\t200\t29\t230\tBacterial small subunit ribosomal RNA\n'
