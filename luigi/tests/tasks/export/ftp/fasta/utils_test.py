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

from Bio.Seq import Seq

from tasks.config import db
from tasks.export.ftp.fasta import utils


def test_can_produce_iterable_from_psql_copy():
    sql = 'select upi, md5 from rna order by id asc limit 5'
    with utils.psql_copy(db(), sql) as raw:
        assert raw.readlines() == [
            'URS0000000001\t6bba097c8c39ed9a0fdf02273ee1c79a\n',
            'URS0000000002\t1fe2f0e3c3a2d6d708ac98e9bfb1d7a8\n',
            'URS0000000003\t7bb11d0572ff6bb42427ce74450ba564\n',
            'URS0000000004\t030c78be0f492872b95219d172e0c658\n',
            'URS0000000005\t030c795b3b5bb84256b0fea3c10ee3aa\n',
        ]


def test_can_turn_tsv_lines_to_seq_records():
    sql = "select upi, md5, 'AAAA' from rna order by id asc limit 3"
    with utils.psql_copy(db(), sql) as tsv_file:
        lines = list(utils.tsv_to_records(tsv_file))
        # Biopython doesn't implement SeqRecord.__eq__, ugh
        assert len(lines) == 3
        assert lines[0].id == 'URS0000000001'
        assert lines[0].description == '6bba097c8c39ed9a0fdf02273ee1c79a'
        assert lines[0].seq == Seq('AAAA')

        assert lines[1].id == 'URS0000000002'
        assert lines[1].description == '1fe2f0e3c3a2d6d708ac98e9bfb1d7a8'
        assert lines[1].seq == Seq('AAAA')

        assert lines[2].id == 'URS0000000003'
        assert lines[2].description == '7bb11d0572ff6bb42427ce74450ba564'
        assert lines[2].seq == Seq('AAAA')


def test_complains_if_no_tsv_lines():
    sql = 'select upi, id, seq_short from rna where id < -1'
    with pytest.raises(ValueError):
        with utils.psql_copy(db(), sql) as out:
            lines = list(utils.tsv_to_records(out))


class SimpleFastaExportBase(utils.FastaExportBase):
    fetch = "select upi, md5, 'AAA' from rna order by id asc limit 5"

    @property
    def filename(self):
        if not hasattr(self, '_handle'):
            self._handle = tempfile.NamedTemporaryFile()
        return self._handle.name


def test_can_produce_correct_file():
    exporter = SimpleFastaExportBase()
    exporter.run()
    with open(exporter.filename, 'r') as out:
        assert out.readlines() == [
            '>URS0000000001 6bba097c8c39ed9a0fdf02273ee1c79a\n',
            'AAA\n',
            '>URS0000000002 1fe2f0e3c3a2d6d708ac98e9bfb1d7a8\n',
            'AAA\n',
            '>URS0000000003 7bb11d0572ff6bb42427ce74450ba564\n',
            'AAA\n',
            '>URS0000000004 030c78be0f492872b95219d172e0c658\n',
            'AAA\n',
            '>URS0000000005 030c795b3b5bb84256b0fea3c10ee3aa\n',
            'AAA\n',
        ]

    exporter._handle.close()
