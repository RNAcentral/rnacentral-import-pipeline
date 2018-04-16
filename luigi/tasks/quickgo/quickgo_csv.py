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

import luigi

from databases.quickgo.parser import parser
from databases.quickgo.data import WRITING_HEADERS

from tasks.config import quickgo

from tasks.utils.fetch import FetchTask


class QuickGoCsv(luigi.Task):
    headers = WRITING_HEADERS

    def requires(self):
        conf = quickgo()
        return [
            FetchTask(remote_path=conf.data_file, local_path=conf.annotations),
        ]

    def output(self):
        return luigi.LocalTarget(quickgo().csv())

    def data(self):
        filename = self.requires()[0].output().fn
        with open(filename, 'w') as raw:
            for go_term in parser(raw):
                yield go_term.as_writeable()
