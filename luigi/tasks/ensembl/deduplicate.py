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
import csv
import shutil
import fileinput
import operator as op
from contextlib import contextmanager

import luigi
from luigi.local_target import atomic_file

from ..config import output
from ..config import ensembl
from ..utils.parameters import CommaGenericFileParameter

from .utils.writers import Output
from .generic import EnsemblSingleFileTask


class DeduplicateOutputType(luigi.Task):  # pylint: disable=R0904
    """
    This is a task that deduplicates a single file type for a single species. For example
    this will deduplicate all Homo sapiens short files.
    """
    filenames = luigi.Parameter()
    output_type = luigi.Parameter()

    @property
    def duplicate_filenames(self):
        """
        Get the list of filenames this task needs to work with.
        """

        files = []
        for task in self.requires():
            target = getattr(task.output(), self.output_type)
            files.append(target.fn)
        return files

    @property
    def unique_columns(self):
        return []

    @contextmanager
    def duplicate_fileobjs(self):
        """
        This is a simple contextmanager to ensure that all file handles are
        closed even if there is an error part way through processing them. It
        will yield a fileinput.input object of all files this task has to
        deduplicate.
        """

        try:
            raw = fileinput.input(files=self.duplicate_filenames)
            yield raw
        finally:
            raw.close()

    def requires(self):
        for filename in self.filenames.split(','):
            yield EnsemblSingleFileTask(input_file=filename)

    def output(self):
        out = Output.build(output().base, 'ensembl', 'dedup')
        final = getattr(out, self.output_type)
        return luigi.LocalTarget(final.fn)

    def copy_files(self, out):
        """
        This will copy the contents of the given files to the output
        filehandle.
        """

        for filename in self.duplicate_filenames:
            with open(filename, 'rb') as raw:
                shutil.copyfileobj(raw, out)

    def dedup_files(self, out, unique_columns):
        """
        This will deduplicate all csv files from self.files() according to the
        columns in unique_columns and write the result to the
        """

        seen = set()
        key = op.itemgetter(*unique_columns)
        writer = csv.writer(out)
        with self.duplicate_fileobjs() as raw:
            reader = csv.reader(raw)
            for row in reader:
                value = key(row)
                if value not in seen:
                    writer.writerow(row)
                    seen.add(value)

    def cleanup(self):
        """
        This will delete all files this has deduplicated.
        """

        for filename in self.duplicate_filenames:
            os.remove(filename)

    def run(self):
        unique_columns = self.unique_columns
        with atomic_file(self.output().fn) as out:
            if not unique_columns:
                self.copy_files(out)
            else:
                self.dedup_files(out, unique_columns)

        # if ensembl().cleanup:
        #     self.cleanup()
