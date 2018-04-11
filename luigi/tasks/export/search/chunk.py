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
from luigi import LocalTarget

from tasks.config import db
from tasks.config import output

from rnacentral.search import exporter


class SearchChunkTask(luigi.Task):  # pylint: disable=R0904
    """
    This is a task that will create an xml export for the given range of UPI's.
    """

    min = luigi.IntParameter()
    max = luigi.IntParameter()

    def output(self):
        config = output()
        filepattern = 'xml4dbdumps__{min}__{max}.xml'.format(
            min=self.min,
            max=self.max,
        )
        filename = os.path.join(config.search_files, filepattern)
        return LocalTarget(filename)

    def run(self):
        with self.output().open('w') as raw:
            results = exporter.range(db(), self.min, self.max)
            exporter.write(raw, results)
