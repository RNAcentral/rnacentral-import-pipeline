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

import gzip
import operator as op

from databases.quickgo import data
from databases.quickgo.parser import parser

from tasks.config import quickgo

from tasks.utils.fetch import FetchTask

from tasks.utils.writers import CsvOutput
from tasks.utils.writers import MultiCsvOutput


class QuickGoData(object):

    def requires(self):
        conf = quickgo()
        return [
            FetchTask(remote_path=conf.data_file, local_path=conf.local_copy),
        ]

    def output(self):
        conf = quickgo()
        return MultiCsvOutput.build(
            annotations=CsvOutput(
                conf.csv,
                data.ANNOTATION_HEADER,
                op.methodcaller('writeable'),
            ),

            publications=CsvOutput(
                conf.publications,
                data.PUB_HEADER,
                op.methodcaller('writeable_publications'),
            ),

            eco=CsvOutput(
                conf.eco_terms,
                data.ECO_HEADER,
                op.methodcaller('writeable_eco_codes'),

            ),

            go_terms=CsvOutput(
                conf.go_terms,
                data.GO_HEADER,
                op.methodcaller('writeable_go_terms'),
            ),
        )

    def run(self):
        filename = self.requires()[0].output()
        with gzip.open(filename, 'r') as raw:
            self.output().populate(parser(raw))
