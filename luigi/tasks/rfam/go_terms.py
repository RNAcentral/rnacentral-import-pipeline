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

import operator as op

import luigi

from databases.rfam import utils

from ontologies import data
from ontologies import helpers as ont

from tasks.config import rfam
from tasks.utils.writers import CsvOutput


class RfamGoTerms(luigi.Task):
    def output(self):
        return CsvOutput(
            rfam().go_terms,
            data.HEADERS,
            op.methodcaller('writeable'),
        )

    def data(self):
        terms = set()
        for family in utils.load_families():
            for (go_term_id, _) in family.go_terms:
                if go_term_id not in terms:
                    terms.add(go_term_id)
                    yield ont.term(go_term_id)

    def run(self):
        self.output().populate(self.data())
