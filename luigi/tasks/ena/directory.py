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
import gzip

import luigi

from databases.ena.parsers import parse_with_mapping_files

from tasks.config import ena
from tasks.config import output
from tasks.utils.fetch import FetchTask
from tasks.utils.entry_writers import Output

from .copy import CopyNcr


class EnaDirectory(luigi.Task):
    input_dir = luigi.Parameter()

    def requires(self):
        yield CopyNcr(ncr=self.input_dir)

        conf = ena()
        for database in ena().tpa_databases:
            local_path = conf.input_tpa_file(database)
            remote_path = conf.raw_tpa_url(database)
            yield FetchTask(remote_path=remote_path, local_path=local_path)

    def output(self):
        prefix = os.path.basename(self.input_dir)
        return Output.build(output().base, 'ena', prefix)

    def handles(self):
        """
        Produce an interable for all compressed non-coding product files that
        this task imports.
        """

        base_dir = CopyNcr(ncr=self.input_dir).output().fn
        files = os.listdir(base_dir)
        for filename in files:
            filename = os.path.join(self.input_dir, filename)
            with gzip.open(filename, 'rb') as raw:
                yield raw

    def run(self):
        files = ena().all_tpa_files()
        with self.output().writer() as writer:
            for handle in self.handles():
                for entry in parse_with_mapping_files(handle, files):
                    if entry.is_valid():
                        writer.write(entry)
