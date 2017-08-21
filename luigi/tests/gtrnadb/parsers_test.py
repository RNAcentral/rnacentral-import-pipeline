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

import unittest as ut

from gtrnadb import parsers


class ExonTest(ut.TestCase):
    def test_it_builds_all_exons(self):
        assert len(self.exons(1)) == 1
        # tRNA-Pro-AGG-1-1 ASM1810v1
        # assert len(self.exons(


class EntryTest(ut.TestCase):
    def test_it_can_generate_all_entries(self):
        assert len(self.entries()) == 1000
