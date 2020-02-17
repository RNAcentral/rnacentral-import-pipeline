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
import operator as op
import itertools as it

from . import fetch

def parse_file(handle):
    reader = csv.reader(handle)
    for line in reader:
        yield fetch.term(line[0])


def process_term_file(term_handle, output):
    data = parse_file(term_handle)
    data = map(op.methodcaller('writeables'), data)
    data = it.chain.from_iterable(data)
    writer = csv.writer(output, quoting=csv.QUOTE_ALL)
    writer.writerows(data)
