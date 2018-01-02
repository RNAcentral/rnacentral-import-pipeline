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
import abc

import luigi
from luigi.local_target import atomic_file

from Bio import SeqIO

from tasks.config import export

from .active import ActiveFastaExport

NHMMER_PATTERN = re.compile('^[ABCDGHKMNRSTVWXYU]+$', re.IGNORECASE)


class NHmmerExportBase(luigi.Task):
    """
    This is the base class that contains most of the logic for writing out an
    nhmmer fasta file. It handles the requirements, output and the actual
    writing. Subclasses need only implement the sequences method.
    """
    __metaclass__ = abc.ABCMeta

    filename = None

    def requires(self):
        return ActiveFastaExport()

    def output(self):
        return luigi.LocalTarget(export().nhmmer(self.filename))

    @abc.abstractmethod
    def is_valid_record(self, record):
        pass

    def sequences(self, filename):
        """
        This should create an interable of the SeqRecords to write out.
        """
        for record in SeqIO.parse(filename, 'fasta'):
            if self.is_valid_record(record):
                yield record

    def run(self):
        filename = self.output().fn
        try:
            os.makedirs(os.path.dirname(filename))
        except:
            pass

        active_file = ActiveFastaExport().output().fn
        with atomic_file(filename) as out:
            SeqIO.write(self.sequences(active_file), out, "fasta")


class NHmmerIncludedExport(NHmmerExportBase):
    """
    This creates a fasta file of all sequences which are included in the nhmmer
    database. This works by parsing the active export files and then producing
    a new one which only contains sequences that contain the accepted charaters
    only.
    """
    filename = 'rnacentral_nhmmer.fasta'

    def is_valid_record(self, record):
        return bool(NHMMER_PATTERN.match(str(record.seq)))


class NHmmerExcludedExport(NHmmerExportBase):
    """
    This creates a fasta file of all sequences which are not included in the
    nhmmer database. This is based off all sequences in the active fasta export
    that contain characters which are not part of the allowed characters.
    """
    filename = 'rnacentral_nhmmer_excluded.fasta'

    def is_valid_record(self, record):
        return not bool(NHMMER_PATTERN.match(str(record.seq)))
