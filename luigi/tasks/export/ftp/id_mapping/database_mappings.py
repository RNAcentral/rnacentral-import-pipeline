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

import csv

import luigi

from tasks.config import export

from .id_mapping import IdMapping


class DatabaseSpecificMappings(luigi.Task):

    def requires(self):
        return IdMapping()

    def output(self):
        pass

    def run(self):
        handles = {}
        with open(IdMapping().output().fn, 'rb') as raw:
            reader = csv.reader(raw, delimiter='\t')
            for row in reader:
                db_name = row[1].lower()
                if db_name not in handles:
                    filename = export().database_mappings(db_name + '.tsv')
                    handle = open(filename)
                    handles[db_name] = (handle, csv.writer(handle))
                handles[db_name][1].writerow(row)

        for (handle, _) in handles.values():
            handle.close()
