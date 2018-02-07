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
from tasks.utils.entry_writers import Output

from . import utils


class SingleEnaFile(luigi.Task):
    """
    Import a Single NCR file from ENA.
    """

    input_file = luigi.Parameter()

    def requires(self):
        yield utils.copy_ncr_task(self.input_file)
        for task in utils.tpa_tasks():
            yield task

    def output(self):
        prefix = os.path.basename(self.input_file)
        return Output.build(output().base, 'ena', prefix)

    def run(self):
        files = ena().all_tpa_files()
        with self.output().writer() as writer:
            with gzip.open(self.input_file, 'rb') as handle:
                writer.write_all(parse_with_mapping_files(handle, files))
