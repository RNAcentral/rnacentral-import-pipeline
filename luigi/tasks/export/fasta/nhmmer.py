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
import re

import luigi
from luigi.local_target import atomic_file

from Bio import SeqIO

from tasks.config import export

from .active import ActiveFastaExport

NHMMER_PATTERN = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


class NHmmerExportBase(luigi.Task):
    def requires(self):
        return ActiveFastaExport()

    def sequences(self):
        return []

    def run(self):
        filename = self.output().fn
        print(filename)
        try:
            os.makedirs(os.path.dirname(filename))
        except:
            pass

        with atomic_file(filename) as out:
            SeqIO.write(self.sequences(), out, "fasta")


class NHmmerDBExport(NHmmerExportBase):
    def output(self):
        return luigi.LocalTarget(export().ftp(
            'sequences',
            '.internal',
            'rnacentral_nhmmer.fasta.gz',
        ))

    def sequences(self):
        filename = ActiveFastaExport().output().fn
        for record in SeqIO.parse(filename, 'fasta'):
            if NHMMER_PATTERN.match(str(record.seq)):
                yield record


class NHmmerExcludedExport(NHmmerExportBase):
    def output(self):
        return luigi.LocalTarget(export().ftp(
            'sequences',
            '.internal',
            'rnacentral_nhmmer_excluded.fasta.gz',
        ))

    def sequences(self):
        filename = ActiveFastaExport().output().fn
        for record in SeqIO.parse(filename, 'fasta'):
            if not NHMMER_PATTERN.match(str(record.seq)):
                yield record
