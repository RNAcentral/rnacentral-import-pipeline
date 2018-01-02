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
    sql = 'select upi, md5 from rna order by id desc limit 5'
    with utils.psql_copy(db(), sql) as raw:
        assert raw.readlines() == [
            'URS0000C8E9EF\tfc30780b063695720cc37a9ce1968b7a\n',
            'URS0000C8E9EE\tf98993c13f10b65603df026553964f5e\n',
            'URS0000C8E9ED\tebad5a2c948b840158b6684319a826e3\n',
            'URS0000C8E9EC\te4a59b51c0426dd5fb73684fb5476537\n',
            'URS0000C8E9EB\tdeec63eaa14fe3571f79a3850545be2a\n',
        ]


def test_can_turn_tsv_lines_to_seq_records():
    sql = "select upi, md5, 'AAAA' from rna order by id desc limit 3"
    with utils.psql_copy(db(), sql) as tsv_file:
        lines = list(utils.tsv_to_records(tsv_file))
        # Biopython doesn't implement SeqRecord.__eq__, ugh
        assert len(lines) == 3
        assert lines[0].id == 'URS0000C8E9EF'
        assert lines[0].description == 'fc30780b063695720cc37a9ce1968b7a'
        assert lines[0].seq == Seq('AAAA')

        assert lines[1].id == 'URS0000C8E9EE'
        assert lines[1].description == 'f98993c13f10b65603df026553964f5e'
        assert lines[1].seq == Seq('AAAA')

        assert lines[2].id == 'URS0000C8E9ED'
        assert lines[2].description == 'ebad5a2c948b840158b6684319a826e3'
        assert lines[2].seq == Seq('AAAA')


def test_complains_if_no_tsv_lines():
    sql = 'select upi, id, seq_short from rna where id < -1'
    with pytest.raises(ValueError):
        with utils.psql_copy(db(), sql) as out:
            lines = list(utils.tsv_to_records(out))
            print(lines)


class SimpleFastaExportBase(utils.FastaExportBase):
    fetch = "select upi, md5, 'AAA' from rna order by id desc limit 5"

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
            '>URS0000C8E9EF fc30780b063695720cc37a9ce1968b7a\n',
            'AAA\n',
            '>URS0000C8E9EE f98993c13f10b65603df026553964f5e\n',
            'AAA\n',
            '>URS0000C8E9ED ebad5a2c948b840158b6684319a826e3\n',
            'AAA\n',
            '>URS0000C8E9EC e4a59b51c0426dd5fb73684fb5476537\n',
            'AAA\n',
            '>URS0000C8E9EB deec63eaa14fe3571f79a3850545be2a\n',
            'AAA\n',
        ]

    exporter._handle.close()
