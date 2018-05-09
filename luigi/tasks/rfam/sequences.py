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

from databases.rfam.parser import parse

from tasks.config import rfam
from tasks.config import output

from tasks.utils.fetch import FetchTask
from tasks.utils.entry_writers import Output


class RfamSequenceFile(luigi.Task):  # pylint: disable=W0232,R0904
    """
    Luigi Task for converting Rfam JSON files into csv files
    that can be loaded into the RNAcentral database.
    """
    input_file = luigi.Parameter()

    def requires(self):
        conf = rfam()
        filename = os.path.basename(self.input_file)
        return FetchTask(
            remote_path=self.input_file,
            local_path=conf.raw('sequences', filename),
        )

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'rfam', prefix)

    def run(self):
        """
        Process json file into RNAcentralEntry objects that can be written to
        the output files using standard methods.
        """
        self.output().populate(parse, self.requires().output())
