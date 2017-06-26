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

import attr
import luigi

from rfam import utils
from rfam.csv_writer import CsvWriter


class FamiliesCSV(CsvWriter):
    headers = [
        'id',
        'name',
        'description',
        'clan',
        'seed_count',
        'full_count',
        'length',
        'domain',
        'is_supressed',
        'rna_type',
    ]

    def data(self):
        for family in utils.load_families():
            data = attr.asdict(family)
            data['is_supressed'] = int(data['is_supressed'])
            yield data


if __name__ == '__main__':
    luigi.run(main_task_cls=FamiliesCSV)
