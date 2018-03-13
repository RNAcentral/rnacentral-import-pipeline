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
from rnacentral.psql import PsqlWrapper


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


def test_can_produce_iterable_from_psql_copy():
    sql = 'select upi, md5 from rna order by id desc limit 5'
    psql = PsqlWrapper(db())
    assert list(psql.copy_to_iterable(sql)) == [
        {'upi': 'URS0000C8E9EF', 'md5': 'fc30780b063695720cc37a9ce1968b7a'},
        {'upi': 'URS0000C8E9EE', 'md5': 'f98993c13f10b65603df026553964f5e'},
        {'upi': 'URS0000C8E9ED', 'md5': 'ebad5a2c948b840158b6684319a826e3'},
        {'upi': 'URS0000C8E9EC', 'md5': 'e4a59b51c0426dd5fb73684fb5476537'},
        {'upi': 'URS0000C8E9EB', 'md5': 'deec63eaa14fe3571f79a3850545be2a'},
    ]
