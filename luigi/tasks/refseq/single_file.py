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

from databases.refseq.parsers import parse

from tasks.config import output
from tasks.config import refseq

from tasks.utils.fetch import FetchTask
from tasks.utils.entry_writers import Output


class RefSeqFile(luigi.Task):
    """
    Parse a single RefSeq file for import.
    """

    input_file = luigi.Parameter()

    def requires(self):
        basename = os.path.basename(self.input_file)
        return FetchTask(
            remote_file=self.input_file,
            local_file=refseq.input_file(basename),
        )

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'refseq', prefix)

    def run(self):
        self.output().populate(parse, self.requires().output())
