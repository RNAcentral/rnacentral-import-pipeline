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

import itertools as it

import luigi

from Bio import SeqIO

from tasks.config import export
from tasks.utils.files import atomic_output

from .active import ActiveFastaExport


class ExampleFasta(luigi.Task):

    def output(self):
        return luigi.LocalTarget(export().sequences('example.txt'))

    def requires(self):
        yield ActiveFastaExport()

    def sequences(self):
        """
        Create an iterable of all sequences to export.
        """

        count = export().search_export_size
        with SeqIO.parse(self.requires().output().fn, 'fasta') as records:
            return list(it.islice(records, count))

    def run(self):
        with atomic_output(self.output()) as out:
            SeqIO.write(self.sequences(), out, 'fasta')
