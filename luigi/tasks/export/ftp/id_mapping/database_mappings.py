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

from tasks.config import export

from rnacentral.export.ftp.id_mapping import split_by_database

from .id_mapping import IdMapping


class DatabaseSpecificMappings(luigi.Task):

    def requires(self):
        return IdMapping()

    def output(self):
        return luigi.LocalTarget(export().database_mappings('ena.tsv'))

    def run(self):
        with self.requires().output().open('r') as raw:
            os.makedirs(export().database_mappings())
            split_by_database(raw, export().database_mappings())
