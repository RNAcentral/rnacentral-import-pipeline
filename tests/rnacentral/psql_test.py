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

import os
import tempfile

import pytest

from rnacentral_pipeline.psql import PsqlWrapper


def test_can_copy_statement_to_iterable():
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    iterable = psql.copy_to_iterable("select upi from rna order by id asc limit 3")
    assert list(iterable) == [
        {'upi': 'URS0000000001'},
        {'upi': 'URS0000000002'},
        {'upi': 'URS0000000003'},
    ]


def test_can_run_command_and_provide_output():
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    with psql.command('explain select upi from rna') as out:
        lines = [line.rstrip('\n') for line in out.readlines()]
        from pprint import pprint
        pprint(lines)
        assert lines == [
            '                                      QUERY PLAN                                       ',
            '---------------------------------------------------------------------------------------',
            ' Index Only Scan using rna_pkey on rna  (cost=0.56..1229364.79 rows=14468335 width=14)',
            '(1 row)',
            '',
        ]


def test_can_produce_iterable_from_psql_copy():
    sql = 'select upi, md5 from rna order by id asc limit 5'
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    assert list(psql.copy_to_iterable(sql)) == [
        {'upi': 'URS0000000001', 'md5': '6bba097c8c39ed9a0fdf02273ee1c79a'},
        {'upi': 'URS0000000002', 'md5': '1fe2f0e3c3a2d6d708ac98e9bfb1d7a8'},
        {'upi': 'URS0000000003', 'md5': '7bb11d0572ff6bb42427ce74450ba564'},
        {'upi': 'URS0000000004', 'md5': '030c78be0f492872b95219d172e0c658'},
        {'upi': 'URS0000000005', 'md5': '030c795b3b5bb84256b0fea3c10ee3aa'},
    ]


def test_it_can_format_a_query():
    sql = 'select upi, md5 from rna order by id asc limit {count}'
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    assert list(psql.copy_to_iterable(sql, count=5)) == [
        {'upi': 'URS0000000001', 'md5': '6bba097c8c39ed9a0fdf02273ee1c79a'},
        {'upi': 'URS0000000002', 'md5': '1fe2f0e3c3a2d6d708ac98e9bfb1d7a8'},
        {'upi': 'URS0000000003', 'md5': '7bb11d0572ff6bb42427ce74450ba564'},
        {'upi': 'URS0000000004', 'md5': '030c78be0f492872b95219d172e0c658'},
        {'upi': 'URS0000000005', 'md5': '030c795b3b5bb84256b0fea3c10ee3aa'},
    ]


def test_can_write_query_to_a_file():
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    with tempfile.NamedTemporaryFile() as tmp:
        psql.write_query(tmp, 'select upi, md5 from rna order by id asc limit 3')
        tmp.flush()
        with open(tmp.name, 'r') as raw:
            assert raw.readlines() == [
                'upi,md5\n',
                'URS0000000001,6bba097c8c39ed9a0fdf02273ee1c79a\n',
                'URS0000000002,1fe2f0e3c3a2d6d708ac98e9bfb1d7a8\n',
                'URS0000000003,7bb11d0572ff6bb42427ce74450ba564\n',
            ]


def test_can_write_a_tsv_file():
    psql = PsqlWrapper(os.environ['PGDATABASE'])
    with tempfile.NamedTemporaryFile() as tmp:
        psql.write_query(tmp, 'select upi, md5 from rna order by id asc limit 3',
                         use='tsv')
        tmp.flush()
        with open(tmp.name, 'r') as raw:
            assert raw.readlines() == [
                'URS0000000001\t6bba097c8c39ed9a0fdf02273ee1c79a\n',
                'URS0000000002\t1fe2f0e3c3a2d6d708ac98e9bfb1d7a8\n',
                'URS0000000003\t7bb11d0572ff6bb42427ce74450ba564\n',
            ]
