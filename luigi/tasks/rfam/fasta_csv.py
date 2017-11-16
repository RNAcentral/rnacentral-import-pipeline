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

from datetime import date

from Bio import SeqIO

from tasks.config import rfam
from tasks.utils.csv_writer import CsvWriter


class RfamFastaCSV(CsvWriter):
    headers = [
        'upi',
        'date',
    ]

    def data(self):
        filename = rfam().fasta
        for record in SeqIO.parse(filename, 'fasta'):
            yield {
                'upi': record.id,
                'date': date.today()
            }
