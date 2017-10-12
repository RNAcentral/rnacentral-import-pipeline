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
import luigi

from tasks.config import output

from .utils.pgload_sequences import PGLoadSequences


class SplitFiles(luigi.Task):
    """
    Split a file into multiple chunks.
    """
    directory = luigi.Parameter()
    extension = luigi.Parameter(default='csv')
    prefix = luigi.Parameter(default='chunk_')

    def requires(self):
        """
        Get the file that needs to be split.
        """
        return MergeFiles(directory=self.directory)

    def run(self):
        """
        Split file into chunks of up to `chunk_size` named
        chunk_00, chunk_01 etc. Lines are preserved.
        """
        cmd = 'cd {path} && split -dC {chunk_size} {merged_file} {prefix}'.format(
            path=os.path.join(output().base, self.directory),
            chunk_size=1024 * 1000 * 1000,
            merged_file=self.input().fn,
            prefix=self.prefix
        )
        os.system(cmd)

    def output(self):
        """
        Check that at least one chunk file exists.
        """
        filename = os.path.join(output().base, self.directory, 'chunk_00')
        return luigi.LocalTarget(filename)


class MergeFiles(luigi.Task):
    """
    Merge all files found in a specified folder into one.
    The taks does not use luigi.ExternalProgramTask because the commands
    are combined using Unix pipes.
    """
    directory = luigi.Parameter()
    extension = luigi.Parameter(default='csv')

    def run(self):
        """
        List all files using `find`.
        Note that `ls` cannot handle folders with an excessive number of files.
        """
        cmd = "find {directory} -type f -name '*.{extension}' | xargs cat > {output}".format(
            directory=os.path.join(output().base, self.directory),
            extension=self.extension,
            output=self.output().fn
        )
        os.system(cmd)

    def output(self):
        """
        Path to the merged file.
        """
        filename = 'merged.%s' % self.extension
        filepath = os.path.join(output().base, self.directory, filename)
        return luigi.LocalTarget(filepath)


class LoadSequences(luigi.WrapperTask):  # pylint: disable=R0904
    """
    This will load sequences into the database. It will load either the short,
    or long sequences as needed. By default, and when type is all, it will load
    botht types of sequences.
    """

    database = luigi.Parameter(default='all')
    type = luigi.ChoiceParameter(
        choices=['short', 'long', 'all'],
        default='all',
    )

    def requires(self):
        if self.type == 'all' or self.type == 'short':
            yield PGLoadSequences(database=self.database, type='short')
        if self.type == 'all' or self.type == 'long':
            yield PGLoadSequences(database=self.database, type='long')
