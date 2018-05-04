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

from tasks.confg import export


class ExportJsonRange(luigi.Task):
    start = luigi.IntParameter()
    stop = luigi.IntParameter()

    def output(self):
        filename = 'ensembl-export-{min}-{max}.json'.format(
            min=self.start,
            max=self.stop,
        )
        return luigi.LocalTarget(export().json(filename))


    def run(self):
        pass
