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
import urlparse

import luigi

from databases.gtrnadb.parsers import parse

from tasks.config import output
from tasks.config import gtrnadb

from tasks.utils import compressed
from tasks.utils.fetch import Fetch
from tasks.utils.entry_writers import Output


class GtRNAdbJsonToCsv(luigi.Task):  # pylint: disable=R0904
    """
    Parse all GtRNAdb JSON files and produce the files needed to import to the
    database.
    """
    url = luigi.Parameter()

    @property
    def filename(self):
        """
        Determine a filename to use for the url. This is the last part of the
        path in the url.
        """

        urlparts = urlparse.urlsplit(self.url)
        parts = urlparts.path.split('/')
        return parts[-1]

    def requires(self):
        filename = gtrnadb().raw(self.filename)
        return Fetch(remote_path=self.url, local_path=filename)

    def output(self):
        prefix = os.path.basename(self.filename)
        return Output.build(output().base, 'gtrnadb', prefix)

    def expanded_files(self):
        """
        This will expand the downloaded, compressed files and provide a list
        the expanded filenames.
        """

        filename = self.requires().output().fn
        dirname = os.path.dirname(filename)
        return compressed.expand(filename, dirname)

    def run(self):
        """
        Create a generator for all entries in all configured GtRNAdb JSON
        files.
        """

        with self.output().writer() as writer:
            for filename in self.expanded_files():
                writer.write_all(parse(filename))
