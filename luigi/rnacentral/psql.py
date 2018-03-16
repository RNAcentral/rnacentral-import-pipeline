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

import re
import csv
import sys
import tempfile
from contextlib import contextmanager
import subprocess as sp


def query_as_copy(sql, use='csv', **kwargs):
    query = str(sql)
    if kwargs:
        query = query.format(**kwargs)
    query = query.replace('\n', ' ')
    query = re.sub('[ ]+', ' ', query)
    options = "DELIMITER AS ',' CSV HEADER"
    if use == 'tsv':
        options = "FORMAT TEXT"
    return "COPY ({query}) TO STDOUT WITH {options}".format(
        query=query,
        options=options,
    )


class PsqlWrapper(object):
    """
    This is a class that wraps some simple psql calls. The idea is that in some
    cases we can't use the direct connection, but instead we should use psql.
    This provides a wrapper around using psql to copy to a file and then
    reading the results from the file.
    """

    def __init__(self, config):
        self.config = config

    @contextmanager
    def command(self, command, handle=None):
        with (handle or tempfile.NamedTemporaryFile()) as out:
            process = sp.Popen(
                [self.config.psql, '-c', command, self.config.pgloader_url()],
                env={'PGPASSWORD': self.config.password},
                stdout=out,
            )

            if process.wait() != 0:
                raise ValueError('Failed to run psql')

            with open(out.name, 'rb') as readable:
                yield readable

    def write_command(self, handle, command):
        with self.command(command, handle=handle):
            pass

    def write_query(self, handle, sql, **kwargs):
        command = query_as_copy(sql, **kwargs)
        self.write_command(command, handle)

    def copy_to_iterable(self, sql, **kwargs):
        """
        This will dump the results of a the query to a TSV file and then create
        a context manage with a file handle of that file. The file is temporary
        and is deleted once the handler exits.
        """

        command = query_as_copy(sql)
        with self.command(command) as out:
            csv.field_size_limit(sys.maxsize)
            for result in csv.DictReader(out):
                yield result
