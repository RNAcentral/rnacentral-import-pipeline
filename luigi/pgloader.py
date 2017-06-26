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
import abc
from tempfile import NamedTemporaryFile

from luigi.contrib.external_program import ExternalProgramTask


class PGLoader(ExternalProgramTask):
    __metaclass__ = abc.ABCMeta

    def program_args(self):
        return ['pgloader', self.__control_filename]

    @abc.abstractmethod
    def control_file(self):
        pass

    def run(self):
        with NamedTemporaryFile() as out:
            out.write(self.control_file())
            self.__control_filename = out.name
            super(PGLoader, self).run()
