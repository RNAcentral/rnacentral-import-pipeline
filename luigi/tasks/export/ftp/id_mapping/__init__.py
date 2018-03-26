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

import luigi

from tasks.utils.compress import GenericCompressTask

from .readme import Readme
from .id_mapping import IdMapping
from .database_mappings import DatabaseSpecificMappings


class Compress(GenericCompressTask):
    def inputs(self):
        yield IdMapping()

    def requires(self):
        for requirement in super(Compress, self).requires():
            yield requirement
        yield DatabaseSpecificMappings()


class IdExport(luigi.WrapperTask):
    def requires(self):
        yield Readme()
        yield DatabaseSpecificMappings()
        yield Compress()
