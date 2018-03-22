# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from tasks.config import db
from tasks.config import genome_mapping

from tasks.export.ftp.fasta.utils import FastaExportBase

from rnacentral.genome_mapping import genome_mapping as gm


class GetFasta(FastaExportBase):
    taxid = luigi.Parameter(default=9606)

    def output(self):
        return luigi.LocalTarget(genome_mapping().rnacentral_fasta('%i.fa' % self.taxid))

    def records(self):
        return gm.export_rnacentral_fasta2(db(), taxid=self.taxid)
