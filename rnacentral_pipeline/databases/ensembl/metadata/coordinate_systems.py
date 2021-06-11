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

import more_itertools as more

from . import databases as db


def is_numeric(entry):
    return isinstance(entry[4], int)


def top_level_only(data):
    parsed = []
    has_name_rank = False
    key = op.itemgetter("name", "coordinate_system", "attrib_value")
    grouped = it.groupby(data, key)
    for (name, sys, assembly), entries in grouped:
        entries = list(entries)
        attribs = {e["attrib_name"]: e["attrib_value"] for e in entries}

        is_ref = attribs.get("non_ref", None) != "1"
        rank = None
        if "karyotype_rank" in attribs:
            try:
                rank = int(attribs["karyotype_rank"])
            except ValueError:
                has_name_rank = True
                rank = attribs["karyotype_rank"]

        # Some chromosomes have name like X_random which we don't want to
        # consider reference, these never have a karyotype_rank so we can
        # exclude them
        is_ref = is_ref and rank is not None
        parsed.append([name, sys, assembly, int(is_ref), rank])

    if not has_name_rank:
        for entry in parsed:
            yield entry

    else:

        def sort_key(items):
            if items[4] is None:
                return ""
            return items[4]

        rest, numeric = more.partition(is_numeric, parsed)
        numeric = list(numeric)
        rest = list(rest)
        ordered = sorted(rest, key=sort_key) + sorted(numeric, key=sort_key)
        index = -1
        for (name, sys, assembly, is_ref, rank) in ordered:
            updated = rank
            if isinstance(rank, (str, int)):
                index += 1
                updated = index

            yield [name, sys, assembly, is_ref, updated]


def fetch(connections, query_handle):
    results = db.run_queries_across_databases(connections, query_handle)
    for (_, rows) in results:
        for entry in top_level_only(rows):
            yield entry


def write(connections, query, output):
    data = fetch(connections, query)
    csv.writer(output).writerows(data)
