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

import luigi

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
    filenames = CommaGenericFileParameter()
    output_type = luigi.Parameter()

    @property
    def files(self):
        """
        Get the list of filenames this task needs to work with.
        """

        return CommaGenericFileParameter.parse(self.filenames)

    def requires(self):
        for filename in self.files:
            yield EnsemblSingleFileTask(input_file=filename)

    def output(self):
        out = Output.build(output().base, 'ensembl', 'dedup')
        final = getattr(out, self.output_type)
        return luigi.LocalTarget(final)

    def copy_files(self, out):
        """
        This will copy the contents of the given files to the output
        filehandle.
        """

        with fileinput.input(files=self.files()) as raw:
            shutil.copyfileobj(raw, out)

    def dedup_files(self, out, unique_columns):
        """
        This will deduplicate all csv files from self.files() according to the
        columns in unique_columns and write the result to the
        """

        seen = set()
        key = op.itemgetter(unique_columns)
        writer = csv.writer(out)
        with fileinput.input(files=self.files()) as raw:
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

        for filename in self.files():
            os.remove(filename)

    def run(self):
        task = EnsemblSingleFileTask(input_file=self.files[0])
        out = task.output()
        output_type = getattr(output, self.output_type)
        klass = output_type.validator.type
        with open(output_type.fn, 'wb') as out:
            unique_columns = getattr(klass, 'unique_columns', None)
            if not unique_columns:
                self.copy_files(out)
            else:
                self.dedup_files(out, unique_columns)

        if ensembl().cleanup:
            self.cleanup()
