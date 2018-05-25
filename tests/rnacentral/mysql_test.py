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

from tasks.config import rfam

from rnacentral.mysql import MysqlWrapper


def test_can_run_a_query():
    mysql = MysqlWrapper(rfam())
    assert list(mysql.query('select * from clan limit 1')) == [{
        'clan_acc': 'CL00001',
        'id': 'tRNA',
        'previous_id': None,
        'description': 'tRNA clan',
        'author': 'Gardner PP',
        'comment': (
            'The tRNA clan contains the RNA families tRNA and tmRNA.'
            'Homology between these families has been established in the'
            'published literature [1-5]'
        ),
        'created': '2014-06-02 15:17:41',
        'updated': '2014-06-02 15:18:18',
    }]


def test_can_write_a_query_to_file():
    mysql = MysqlWrapper(rfam())
    with tempfile.NamedTemporaryFile() as tmp:
        mysql.write_query(tmp, 'select clan_acc, id from clan limit 3')
        tmp.flush()
        with open(tmp.name, 'r') as raw:
            assert raw.readlines() == [
                'clan_acc\tid',
                'CL00001\ttRNA',
                'CL00002\tRNaseP',
                'CL00003\tSRP',
            ]


def test_can_run_a_command():
    mysql = MysqlWrapper(rfam())
    with mysql.command('expalin select id from clan') as out:
        lines = [l.rstrip('\n') for l in out.readlines()]
        assert lines == [
            'id\tselect_type\ttable\ttype\tpossible_keys\tkey\tkey_len\tref\trows\tExtra',
            '1\tSIMPLE\tclan\tALL\tNULL\tNULL\tNULL\tNULL\t112\tNULL',
        ]
