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

from ensembl import data
from ensembl.helpers import bio as helpers


class Base(ut.TestCase):  # pylint: disable=R0904
    filename = None
    importer_class = None

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
            elif not helpers.is_ncrna(feature):
                continue

            if helpers.is_transcript(feature):
                key = helpers.transcript(feature)
            if not key:
                continue
            cls.features[key] = feature

    def setUp(self):
        self.importer = None
        if self.importer_class and self.filename:
            self.importer = self.importer_class()

    def data(self):
        with open(self.filename, 'rb') as raw:
            for entry in self.importer.data(raw):
                yield entry

    def summary_of(self, key):
        summary = data.Summary(sequence=self.record.seq)
        return summary.update_gene_info(self.features[key])

    def entries_for(self, feature_key):
        feature = self.features[feature_key]
        summary = self.summary_of(helpers.gene(feature))
        entries = self.importer.rnacentral_entries(self.record, summary,
                                                   feature)
        return list(entries)

    def entry_for(self, feature_key):
        entries = self.entries_for(feature_key)
        assert len(entries) == 1
        return entries[0]

    def entries_from(self, database, feature_key):
        entries = self.entries_for(feature_key)
        return [e for e in entries if e.database == database]

    def entry_from(self, database, feature_key):
        entries = self.entries_from(database, feature_key)
        assert len(entries) == 1
        return entries[0]
