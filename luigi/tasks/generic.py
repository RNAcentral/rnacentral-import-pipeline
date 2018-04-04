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

from databases.generic.parser import parse

from tasks.config import output
from tasks.config import generic
from tasks.utils.fetch import FetchTask
from tasks.utils.entry_writers import Output


class GenericDatabase(luigi.Task):  # pylint: disable=too-many-public-methods
    """
    This tasks handles importing data for a database that uses our generic JSON
    schema based import. It will fetch the data and then create CSV's for
    pgloader from it.
    """

    input_file = luigi.Parameter()

    def requires(self):
        config = generic()
        return FetchTask(
            remote_path=self.input_file,
            local_path=config.raw(self.input_file),
        )

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'generic', prefix)

    def run(self):
        self.output().populate(parse, self.requires().output())
