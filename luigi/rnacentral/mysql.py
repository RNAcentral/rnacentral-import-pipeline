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
        """
        Runs a command and yields a filehandle that contains the output.
        """

        with tempfile.NamedTemporaryFile() as out:
            cmd = self.mysql[
                '--host', self.config.host,
                '--port', str(self.config.port),
                '--user', self.config.user,
                '--database', self.config.db_name,
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
        """
        This will write the results of a query to the given filehandle.
        """
        self.write_command(handle, query)

    def query(self, query):
        """
        Run a query and return an iterable of the results.
        """

        with self.command(query) as out:
            csv.field_size_limit(sys.maxsize)
            for result in csv.DictReader(out):
                yield result
