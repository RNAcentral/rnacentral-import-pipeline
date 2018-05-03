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

import csv
import sys
import shutil
import tempfile
from contextlib import contextmanager

from plumbum import local


class MysqlWrapper(object):
    """
    A wrapper around the command line mysql program.
    """

    def __init__(self, config):
        self.config = config
        self.mysql = local['mysql']

    @contextmanager
    def command(self, command):
        with tempfile.NamedTemporaryFile() as out:
            cmd = self.mysql[
                '--host', self.config.hostname,
                '--port', str(self.config.port),
                '--user', self.config.username,
                '--database', self.config.database,
            ]

            ((cmd << command) > out)()
            out.flush()

            with open(out.name, 'rb') as readable:
                yield readable

    def write_command(self, handle, command):
        """
        Write the command to the given file handle.
        """

        with self.command(command) as tmp:
            shutil.copyfileobj(tmp, handle)

    def write_query(self, handle, query):
        self.write_command(handle, query)

    def write_table(self, handle, table):
        """
        Download the specific table using MySQL to a specific file. The input URI
        should be of the form:
            mysql://rfamro@mysql-rfam-public.ebi.ac.uk:4497/Rfam?database_link

        that is:
            mysql://username@host:port/database?table
        """
        return self.write_query(handle, 'select * from %s' % table)

    def query(self, query):
        """
        Run a query and return an iterable of the results.
        """

        with self.command(query) as out:
            csv.field_size_limit(sys.maxsize)
            for result in csv.DictReader(out):
                yield result
