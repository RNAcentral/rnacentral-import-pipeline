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

from tasks.config import mirbase
from tasks.generic import GenericDatabase


class MirBase(luigi.WrapperTask):
    """
    Import all Lncipedia data. The configuration file must include the path to
    the JSON file to import. The path may be a url to download from as well,
    but it may not be compressed.
    """

    def requires(self):
        yield GenericDatabase(input_file=mirbase().json_file)
