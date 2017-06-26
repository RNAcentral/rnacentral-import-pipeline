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

from rfam.pgload_hits import PGLoadHits
from rfam.pgload_clans import PGLoadClans
from rfam.pgload_families import PGLoadFamilies


class Rfam(luigi.Task):
    def output(self):
        """
        Create an iterator of all outputs that this will produce.

        Yields
        ------
        An Output object for each output this will produce.
        """
        for requirement in self.requires():
            yield requirement.output()

    def requires(self):
        yield PGLoadClans()
        yield PGLoadFamilies()
        yield PGLoadHits()


if __name__ == '__main__':
    luigi.run(main_task_cls=Rfam)
