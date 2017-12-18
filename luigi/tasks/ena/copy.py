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
from glob import glob

import luigi
from luigi.contrib.external_program import ExternalProgramTask

from tasks.config import ena


class CopyNcr(ExternalProgramTask):
    ncr = luigi.Parameter()

    def output(self):
        prefix = os.path.commonprefix([self.ncr, ena().base])
        if not os.path.isdir(prefix):
            raise ValueError("Could not find common path")
        basename = os.path.relpath(self.ncr, prefix)
        return luigi.LocalTarget(ena().input_ncr_file(basename))

    def copy_file(self):
        try:
            os.makedirs(os.path.dirname(self.output().fn))
        except:
            pass
            return [
                'cp',
                self.ncr,
                self.output().fn,
            ]

    def copy_directory(self):
        try:
            os.makedirs(self.output().fn)
        except:
            pass

        cmd = ['cp', '-r']
        filenames = glob(os.path.join(self.ncr, '*.ncr.gz'))
        if not filenames:
            return []

        cmd.extend(filenames)
        cmd.append(self.output().fn)
        return cmd

    def program_args(self):
        if os.path.isfile(self.ncr):
            return self.copy_file()

        if os.path.isdir(self.ncr):
            return self.copy_directory()

        raise ValueError("Must give file or directory: " + self.ncr)
