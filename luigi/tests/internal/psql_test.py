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

from tasks.config import db
from internal.psql import PsqlWrapper


def test_can_copy_statement_to_iterable():
    psql = PsqlWrapper(db())
    iterable = psql.copy_to_iterable("select upi from rna order by id limit 3")
    assert list(iterable) == [
        {'upi': 'URS0000000001'},
        {'upi': 'URS0000000002'},
        {'upi': 'URS0000000003'},
    ]


def test_can_run_command_and_provide_output():
    psql = PsqlWrapper(db())
    with psql.command('explain select upi from rna') as out:
        lines = [line.rstrip('\n') for line in out.readlines()]
        assert lines == [
            '                           QUERY PLAN                            ',
            '-----------------------------------------------------------------',
            ' Seq Scan on rna  (cost=0.00..1258736.79 rows=13167479 width=14)',
            '(1 row)',
            '',
        ]
