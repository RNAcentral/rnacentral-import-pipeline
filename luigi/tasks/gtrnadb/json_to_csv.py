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

import luigi

from gtrnadb.parsers import parse

from tasks.config import output

from tasks.utils.parameters import PathParameter
from tasks.utils.entry_writers import Output


class GtRNAdbJsonToCsv(luigi.Task):  # pylint: disable=R0904
    """
    Parse all GtRNAdb JSON files and produce the files needed to import to the
    database.
    """
    input_file = PathParameter()

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'gtrnadb', prefix)

    def run(self):
        """
        Create a generator for all entries in all configured GtRNAdb JSON files.
        """

        with self.output().writer() as writer:
            for entry in parse(self.input_file):
                if entry.is_valid():
                    writer.write(entry)
