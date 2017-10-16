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
# import fileinput
import subprocess as sp
# from contextlib import contextmanager

import luigi
from luigi.local_target import atomic_file
# from luigi.contrib.external_program import ExternalProgramTask

from ..config import output

from .generic import EnsemblFile

from tasks.utils.entry_writers import Output


class DeduplicateOutputType(luigi.Task):  # pylint: disable=R0904
    """
    This is a task that deduplicates a single file type for a single species. For example
    this will deduplicate all Homo sapiens short files.
    """
    species_name = luigi.Parameter()
    filenames = luigi.Parameter()
    output_type = luigi.Parameter()

    @property
    def unique_columns(self):
        """"
        This will look at the writer for the given type and call the method to
        determine what columns in the output should be used to detect the
        unique rows.
        """

        out = Output.build(output().base, 'ensembl', self.species_name + '-dedup')
        final = getattr(out.writer(), self.output_type)
        return final.unique_columns()

    def requires(self):
        for filename in self.filenames.split(','):
            yield EnsemblFile(input_file=filename)

    def output(self):
        out = Output.build(output().base, 'ensembl', self.species_name + '-dedup')
        final = getattr(out, self.output_type)
        return luigi.LocalTarget(final.fn)

    def filename_pattern(self):
        return os.path.join(
            output().duplicated_ensembl_path(),
            self.output_type,
            'ensembl*{species}*.csv'.format(species=self.species_name),
        )

    def run(self):
        """
        This will deduplicate all csv files from self.files() according to the
        columns in unique_columns and write the result to the
        """

        start = None
        stop = None
        unique_columns = self.unique_columns
        if len(unique_columns) == 1:
            start = unique_columns[0] + 1
            stop = unique_columns[0] + 1
        elif len(unique_columns) == 2:
            start, stop = unique_columns
            start += 1
            stop += 1
        else:
            raise ValueError("Unique columns must be a tuple of start, stop")

        command = "sort -u -t, -k{start},{stop} {pattern}".format(
            start=start,
            stop=stop,
            pattern=self.filename_pattern(),
        )

        with atomic_file(self.output().fn) as out:
            sp.check_call(command, stdout=out, shell=True)
