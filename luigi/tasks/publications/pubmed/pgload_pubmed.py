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

from tasks.utils.pgloader import PGLoader

from .publication_csv import PubmedData

CONTROL_FILE = """
"""


class PubmedLoader(PGLoader):
    def requires(self):
        return [
            PubmedData(),
        ]

    def control_file(self):
        filename = self.requires()[0].output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='ref_pubmed'),
            search_path=self.db_search_path(),
        )
