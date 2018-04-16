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

from ontologies.eco import to_load
from ontologies.eco import TermSources

from tasks.config import quickgo

from tasks.utils.fetch import FetchTask
from tasks.utils.csv_writer import CsvWriter


class EcoCodeCSV(CsvWriter):
    headers = [
        'eco_term_id',
        'name',
        'description',
    ]

    def requires(self):
        conf = quickgo()
        return [
            FetchTask(remote_path=conf.data_file, local_path=conf.annotations),
        ]

    def data(self):
        source = TermSources(
            quickgo_file=self.requires()[0].output().fn,
        )

        for term in to_load(source):
            yield {
                'eco_term_id': term.ontology_id,
                'name': term.name,
                'description': term.definition,
            }
