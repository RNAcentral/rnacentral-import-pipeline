# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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


def only_or_fail(possible):
    assert len(possible) == 1
    return possible.pop()


def parse(handle):
    reader = csv.reader(handle, delimiter='\t')
    grouped = it.groupby(reader, op.itemgetter(0))
    for gene_id, entries in grouped:
        entries = list(entries)
        description = only_or_fail({e[1] for e in entries}).strip()
        symbol = only_or_fail({e[2] for e in entries})
        synonyms = set()
        for entry in entries:
            value = entry[3].replace('"', '')
            if value:
                synonyms.add(value)

        synonyms = ','.join('"%s"' % s for s in synonyms)
        synonyms = '{%s}' % synonyms
        yield [
            'ENSEMBL:%s' % gene_id,
            description,
            symbol,
            synonyms
        ]


def from_file(handle, output):
    writer = csv.writer(output)
    writer.writerows(parse(handle))
