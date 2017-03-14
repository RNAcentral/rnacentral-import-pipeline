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


class Base(ut.TestCase):
    filename = None
    importer_class = None

    @classmethod
    def setUpClass(cls):
        cls.features = {}
        cls.record = SeqIO.read(cls.filename, 'embl')
        for feature in cls.record.features:
            gene = feature.qualifiers.get('gene', [None])[0]
            key = (feature.type, gene)
            cls.features[key] = feature

    def setUp(self):
        self.importer = self.importer_class(self.filename)
