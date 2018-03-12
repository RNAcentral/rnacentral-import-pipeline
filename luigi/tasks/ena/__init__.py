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

from tasks.config import ena

from .single_file import SingleEnaFile
from .directory import EnaDirectory
from .tpa import FetchTPA


class Ena(luigi.WrapperTask):

    def requires(self):
        for database in ena().tpa_databases:
            yield FetchTPA(database=database)

        for entry in ena().raw_ncr_files():
            if os.path.isdir(entry):
                yield EnaDirectory(input_dir=entry)
            elif os.path.isfile(entry):
                yield SingleEnaFile(input_file=entry)
