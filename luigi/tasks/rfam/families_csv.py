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

from databases.rfam import utils
from tasks.utils.csv_writer import CsvWriter


class RfamFamiliesCSV(CsvWriter):
    headers = [
        'id',
        'short_name',
        'long_name',
        'description',
        'clan_id',
        'seed_count',
        'full_count',
        'length',
        'domain',
        'is_suppressed',
        'rna_type',
        'rfam_rna_type',
    ]

    def data(self):
        for family in utils.load_families():
            data = attr.asdict(family)
            data['short_name'] = family.name
            data['long_name'] = family.pretty_name
            data['is_suppressed'] = int(family.is_suppressed)
            data['rfam_rna_type'] = data['rna_type']
            data['rna_type'] = family.guess_insdc()
            yield data
