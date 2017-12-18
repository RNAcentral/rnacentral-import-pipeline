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

from tasks.config import ena
from tasks.utils.http import download

class FetchTPA(luigi.Task):
    database = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(ena().input_tpa_file(self.database))

    def run(self):
        download(ena().raw_tpa_url(self.database), self.output().fn)
