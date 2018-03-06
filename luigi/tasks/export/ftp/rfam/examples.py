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

import subprocess as sp

import luigi
from luigi.local_target import atomic_file

from tasks.config import export

from .rfam_annotations import RfamAnnotations


class RfamAnnotationsExample(luigi.Task):
    def requires(self):
        return RfamAnnotations()

    def output(self):
        filename = export().rfam('example.txt')
        return luigi.LocalTarget(filename)

    def program_args(self):
        return [
            'head',
            '-n',
            str(export().rfam_example_size),
            self.requires().output().fn,
        ]

    def run(self):
        with atomic_file(self.output().fn) as out:
            sp.Popen(self.program_args(), stdout=out)
