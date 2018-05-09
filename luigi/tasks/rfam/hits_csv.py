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

from glob import iglob

import attr

from databases.rfam.infernal_results import parse
from tasks.config import rfam
from tasks.utils.csv_writer import CsvWriter


class RfamHitsCSV(CsvWriter):  # pylint: disable=R0904
    """
    This class will build a CSV file that can be loaded by pgloader from the
    Rfam hits. It uses the rfam_hits configuration value in the files section.
    """

    headers = [
        'seq_name',
        'seq_from',
        'seq_to',
        'strand',
        'rfam_acc',
        'mdl_from',
        'mdl_to',
        'overlap',
        'e_value',
        'score',
    ]

    def data(self):
        for filename in iglob(rfam().hits):
            for hit in parse(filename, clan_competition=True):
                if hit.overlap == 'unique' or hit.overlap == 'best':
                    yield attr.asdict(hit)
