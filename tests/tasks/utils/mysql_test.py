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

from tasks.utils import mysql


def test_can_fetch_a_table():
    with tempfile.NamedTemporaryFile() as tmp:
        mysql.download(
            'mysql://rfamro@mysql-rfam-public.ebi.ac.uk:4497/Rfam?database_link',
            tmp.name,
        )

        with open(tmp.name, 'r') as raw:
            data = raw.readlines()[0:2]
            assert data == [
                'rfam_acc\tdb_id\tcomment\tdb_link\tother_params\n',
                'RF00014\tSO\t\t0000378\tDsrA_RNA\n',
            ]
