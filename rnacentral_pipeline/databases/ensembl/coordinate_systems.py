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


def top_level_only(data):
    grouped = it.groupby(data, op.itemgetter(0, 1, 2))
    for (name, sys, assembly), entries in grouped:
        entries = list(entries)
        attribs = {e[3]: e[4] for e in entries}

        is_ref = attribs.get('non_ref', None) != '1'
        rank = None
        if 'karyotype_rank' in attribs:
            rank = int(attribs['karyotype_rank'])

        # Some chromosomes have name like X_random which we don't want to
        # consider reference, these never have a karyotype_rank so we can
        # exclude them
        is_ref = is_ref and rank is not None
        yield [name, sys, assembly, int(is_ref), rank]


def from_file(handle, output):
    reader = csv.reader(handle, delimiter='\t')
    writer = csv.writer(output)
    writer.writerows(top_level_only(list(reader)))
