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

from databases.rfam import utils
from tasks.utils.csv_writer import CsvWriter


class GoTermsCSV(CsvWriter):
    headers = [
        'go_term_id',
        'name',
    ]

    def data(self):
        seen = set()
        for family in utils.load_families():
            for (go_term_id, name) in family.go_terms:
                if go_term_id not in seen:
                    seen.add(go_term_id)
                    yield {
                        'go_term_id': go_term_id,
                        'name': name,
                    }
