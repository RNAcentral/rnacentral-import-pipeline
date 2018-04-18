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

import tempfile

from tasks.utils import writers

def test_can_create_usable_csv_writer():

    def getter(entry):
        yield {'a': entry['a'] * -1, 'b': entry['b'] * -1}

    with tempfile.NamedTemporaryFile() as tmp:
        writer = writers.CsvOutput(tmp.name, ["a", "b"], getter)
        writer.populate([
            {'a': 1, 'b': 3, 'c': 0},
            {'a': 2, 'b': 4, 'c': 10},
        ])
        tmp.flush()
        with open(tmp.name, 'r') as result:
            assert result.readlines() == [
                '"-1","-3"\n',
                '"-2","-4"\n',
            ]
