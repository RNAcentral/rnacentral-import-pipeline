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

from ontologies import data
from databases.rfam import cross_references as cr

from tasks.config import rfam

from tasks.utils.fetch import FetchTask
from tasks.utils.writers import CsvOutput


class RfamOntologyTerms(luigi.Task):  # pylint: disable=too-many-public-methods
    """
    A task that will parse ontology data from Rfam to produce data to load.
    """

    def requires(self):
        conf = rfam()
        return FetchTask(
            remote_path=conf.remote_table('database_link'),
            local_path=conf.raw('database_link.tsv'),
        )

    def output(self):
        return CsvOutput(
            rfam().go_terms,
            data.HEADERS,
            op.methodcaller('writeables'),
        )

    def run(self):
        with self.requires().output().open('r') as raw:
            self.output().populate(cr.ontology_terms(raw))
