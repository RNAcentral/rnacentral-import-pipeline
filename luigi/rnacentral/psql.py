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
import shutil
from contextlib import contextmanager

from plumbum import local


def query_as_copy(sql, use='csv', **kwargs):
    """
    Turn a sql query string into a a COPY command for psql. This will format
    the query with the given options. If given the keyword argument use='tsv'
    this will produce a tsv file.
    """

    query = str(sql)
    if kwargs:
        query = query.format(**kwargs)
    query = query.replace('\n', ' ')
    query = re.sub('[ ]+', ' ', query)
    options = "WITH DELIMITER AS ',' CSV HEADER"
    if use == 'tsv':
        options = ""
    return "COPY ({query}) TO STDOUT {options}".format(
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
        self.psql = local['psql']

    @contextmanager
    def command(self, command):
        """
        Execute a command and yield a filehandle that stores the results.
        """

        with local.env(PGPASSWORD=self.config.password):
            cmd = self.psql['-c', command, self.config.pgloader_url()]
            obj = cmd.popen()
            yield obj.stdout

    def write_command(self, handle, command):
        """
        This will write the command to the given file handle.
        """

        with self.command(command) as tmp:
            shutil.copyfileobj(tmp, handle)

    def write_query(self, handle, sql, **kwargs):
        """
        This will write the sql query to the given file handle.
        """

        command = query_as_copy(sql, **kwargs)
        self.write_command(handle, command)

    def copy_to_iterable(self, sql, **kwargs):
        """
        This will dump the results of a the query to a TSV file and then create
        a context manage with a file handle of that file. The file is temporary
        and is deleted once the handler exits.
        """

        command = query_as_copy(sql, **kwargs)
        with self.command(command) as out:
            csv.field_size_limit(sys.maxsize)
            for result in csv.DictReader(out):
                yield result
