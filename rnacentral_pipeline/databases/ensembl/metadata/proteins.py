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

import re
import csv
import operator as op
import itertools as it

from . import databases as db


def only_or_fail(possible):
    assert len(possible) == 1
    return possible.pop()


def short_description(entries):
    description = {e['description'] for e in entries}
    description.discard(None)
    if not description:
        return None
    description = only_or_fail(description).strip()
    if description == 'NULL':
        return None
    description = re.sub(r'\s*\[.+\]$', '', description)
    return description


def parse(data):
    grouped = it.groupby(data, op.itemgetter('stable_id'))
    for gene_id, entries in grouped:
        entries = list(entries)
        symbol = only_or_fail({e['display_label'] for e in entries})
        description = short_description(entries)
        synonyms = set()
        for entry in entries:
            value = entry['synonym']
            if value and value != 'NULL':
                synonyms.add(value.replace('"', ''))

        synonyms = ','.join('"%s"' % s for s in synonyms)
        synonyms = '{%s}' % synonyms
        yield [
            'ENSEMBL:%s' % gene_id,
            description,
            symbol,
            synonyms
        ]


def fetch(connections, query_handle):
    results = db.run_queries_across_databases(connections, query_handle)
    for (_, rows) in results:
        for protein in parse(rows):
            yield protein


def write(connections, query, output):
    data = fetch(connections, query)
    csv.writer(output).writerows(data)
