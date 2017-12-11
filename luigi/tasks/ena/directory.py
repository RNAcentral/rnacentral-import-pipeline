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

from databases.ena.parsers import parse

from tasks.config import output
from tasks.utils.entry_writers import Output


class EnaDirectory(luigi.Task):
    input_dir = luigi.Parameter()

    def output(self):
        prefix = os.path.basename(self.input_dir)
        return Output.build(output().base, 'ena', prefix)

    def files(self):
        files = os.listdir(self.input_dir)
        for filename in files:
            if filename == 'fasta':
                continue
            filename = os.path.join(self.input_dir, filename)
            with gzip.open(filename, 'rb') as raw:
                yield raw

    def run(self):
        with self.output().writer() as writer:
            for handle in self.files():
                for entry in parse(handle):
                    if entry.is_valid():
                        writer.write(entry)
