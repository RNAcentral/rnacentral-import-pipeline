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

import psycopg2

import luigi
from luigi.local_target import atomic_file

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tasks.release.utils.db import get_db_connection
from tasks.config import db


class FastaExportBase(luigi.Task):
    table = None
    populate = None
    fetch = None

    def sequences(self, cursor):
        print(self.fetch)
        cursor.execute(self.fetch)
        print('executed')
        for result in cursor:
            print(result)
            yield SeqRecord(
                Seq(result['sequence']),
                id=result['id'],
                description=result['description'],
            )

    def run(self):
        connection = get_db_connection(db(), connect_timeout=20 * 60)
        # populate_table_by_partition(connection, self.table, self.populate)
        # unique = md5(self.table).\
        #     update(self.populate).\
        #     update(self.fetch).\
        #     hexdigest()

        cursor = connection.cursor(
            cursor_factory=psycopg2.extras.DictCursor,
            name=self.__class__.__name__,
        )
        filename = self.output().fn
        try:
            os.makedirs(os.path.basename(filename))
        except:
            pass

        with atomic_file(filename) as out:
            SeqIO.write(self.sequences(cursor), out, "fasta")
        cursor.close()
        connection.close()


def populate_table_by_partition(connection, table, populate):
    """
    This will create and populate some temporary table using the given table
    create statement and populate statement. The population is done by
    partitioned xref table, and the statement will be formatted with dbid prior
    to being run.
    """

    cursor = connection.cursor()
    cursor.execute(table)
    cursor.execute('select id from rnc_database order by id')
    dbids = [r[0] for r in cursor.fetchall()]
    for dbid in dbids:
        cursor.execute(populate.format(dbid=dbid))
