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
from contextlib import contextmanager
import subprocess as sp

from distutils.spawn import find_executable

import luigi
from luigi.local_target import atomic_file

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tasks.config import db
from tasks.config import export


def tsv_to_records(tsv_lines):
    """
    Given a file of TSV entries (id description sequence) produce a generator
    of Bio.SeqRecord.SeqRecord items.
    """

    seen = False
    for line in tsv_lines:
        sid, description, sequence = line.strip().split('\t')
        seen = True
        yield SeqRecord(
            Seq(sequence),
            id=sid,
            description=description,
        )

    if not seen:
        raise ValueError("No lines found")


@contextmanager
def psql_copy(config, sql, psql=None):
    """
    This will dump the results of a the query to a TSV file and then create
    a context manage with a file handle of that file. The file is temporary
    and is deleted once the handler exits.
    """

    if psql is None:
        psql = find_executable('psql')

    with tempfile.NamedTemporaryFile() as out:
        query = sql.replace('\n', ' ')
        command = 'COPY ({query}) TO STDOUT'.format(query=query)
        process = sp.Popen(
            [psql, '-c', command, config.pgloader_url()],
            env={'PGPASSWORD': config.password},
            stdout=out,
        )

        if process.wait() != 0:
            raise ValueError('Failed to run psql')

        with open(out.name, 'rb') as readable:
            yield readable


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


class FastaExportBase(luigi.Task):
    filename = None
    fetch = None

    def output(self):
        return luigi.LocalTarget(export().sequences(self.filename))

    def run(self):
        filename = self.output().fn
        try:
            os.makedirs(os.path.dirname(filename))
        except:
            pass

        with psql_copy(db(), self.fetch) as tsv_file, atomic_file(filename) as out:
            SeqIO.write(tsv_to_records(tsv_file), out, "fasta")
