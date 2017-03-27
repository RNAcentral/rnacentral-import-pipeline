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

import unittest as ut

from Bio import SeqIO
import luigi

from ensembl import helpers


class Base(ut.TestCase):
    filename = None
    importer_class = None

    def target(self):
        return luigi.LocalTarget(path=self.filename)

    def setUp(self):
        self.importer = None
        if self.importer_class and self.filename:
            self.importer = self.importer_class(self.filename)

    def data(self):
        return self.importer.data(self.target())

    def summary_of(self, key):
        summary = self.importer.summary(self.record)
        return self.importer.update_gene_info(summary, self.features[key])

    @classmethod
    def setUpClass(cls):
        cls.features = {}
        if not cls.filename:
            return

        cls.record = SeqIO.read(cls.filename, 'embl')
        for feature in cls.record.features:
            key = None
            if helpers.is_gene(feature):
                key = helpers.gene(feature)
            elif helpers.is_transcript(feature):
                key = helpers.transcript(feature)
            if not key:
                continue
            cls.features[key] = feature
