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
import itertools as it

import pytest
from Bio.Seq import Seq

from rnacentral.export.ftp import fasta
from rnacentral.psql import PsqlWrapper

from tasks.config import db
from tasks.export.ftp.fasta import utils


class SimpleFastaExportBase(utils.FastaExportBase):

    @property
    def filename(self):
        if not hasattr(self, '_handle'):
            self._handle = tempfile.NamedTemporaryFile()
        return self._handle.name

    def records(self):
        psql = PsqlWrapper(db())
        fetch = """
            select
              upi as id,
              md5 as description,
              'AAA' as sequence
            from rna
            order by id asc
            limit 5
        """
        return it.imap(fasta.as_record, psql.copy_to_iterable(fetch))


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
