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

import abc

import luigi

from Bio import SeqIO

from tasks.config import export
from tasks.utils.files import atomic_output


class FastaExportBase(luigi.Task):
    """
    A base class for tasks that export fasta files.
    """
    __metaclass__ = abc.ABCMeta

    filename = None

    def output(self):
        return luigi.LocalTarget(
            export().sequences(self.filename),
            format=luigi.format.Gzip,
        )

    @abc.abstractmethod
    def records(self):
        """
        This should generate an iterable of Bio.SeqRecords to write out.
        """
        return []

    def run(self):
        with atomic_output(self.output()) as out:
            SeqIO.write(self.records(), out, "fasta")
