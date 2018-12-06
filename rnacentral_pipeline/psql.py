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
import json
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


def json_handler(out):
    for row in out:
        row = row.replace('\\\\', '\\')
        yield json.loads(row)


def csv_handler(out):
    csv.field_size_limit(sys.maxsize)
    for result in csv.DictReader(out):
        yield result


class PsqlWrapper(object):
    """
    This is a class that wraps some simple psql calls. The idea is that in some
    cases we can't use the direct connection, but instead we should use psql.
    This provides a wrapper around using psql to copy to a file and then
    reading the results from the file.
    """

    def __init__(self, url):
        self.pgloader_url = url
        self.psql = local['psql']

    @contextmanager
    def command(self, command, options=None):
        """
        Execute a command and yield a filehandle that stores the results.
        """

        cmd = self.psql['-c', command, self.pgloader_url]
        obj = cmd.popen()
        yield obj.stdout

        if obj.returncode:
            raise ValueError(obj.stderr.read())

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
        self.write_command(handle, query_as_copy(sql, **kwargs))

    def copy_file_to_handle(self, filename, handle):
        options = ['-f', filename, self.pgloader_url]
        cmd = self.psql.bound_command(*options)
        obj = cmd.popen(stdout=handle)
        if obj.returncode:
            raise ValueError(obj.stderr.read())
        obj.wait()

    def copy_file_to_iterable(self, filename, json=False, **kwargs):
        handler = csv_handler
        if json:
            handler = json_handler

        options = ['-f', filename, self.pgloader_url]
        cmd = self.psql.bound_command(*options)
        obj = cmd.popen()
        for result in handler(obj.stdout):
            yield result

        if obj.returncode:
            raise ValueError(obj.stderr.read())

    def copy_to_iterable(self, sql, json=False, **kwargs):
        """
        This will dump the results of a the query to a TSV file and then create
        a context manage with a file handle of that file. The file is temporary
        and is deleted once the handler exits.
        """

        handler = csv_handler
        if json:
            kwargs['tsv'] = True
            handler = json_handler

        query = query_as_copy(sql, **kwargs)
        with self.command(query) as out:
            return handler(out)
