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

from tasks.config import rgd
from tasks.utils.fetch import fetch

from databases.rgd import helpers


class RgdFetchGff(luigi.Task):
    organism = luigi.Parameter()

    @property
    def remote_files(self):
        if not hasattr('__gff_files'):
            self.__gff_files = helpers.gff_files(rgd().host, self.organism)
        return self.__gff_files

    def output(self):
        basename = os.path.basename(self.remote_files)
        filename = rgd().raw(self.organism, basename)
        return luigi.LocalTarget(filename)

    def run(self):
        for filename in self.remote_files:
            fetch(filename, self.output().fn)
