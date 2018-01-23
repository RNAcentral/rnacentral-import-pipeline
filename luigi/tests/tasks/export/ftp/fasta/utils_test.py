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

from internal.export.ftp import fasta
from internal.psql import PsqlWrapper

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
            order by id desc
            limit 5
        """
        return it.imap(fasta.as_record, psql.copy_to_iterable(fetch))


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
