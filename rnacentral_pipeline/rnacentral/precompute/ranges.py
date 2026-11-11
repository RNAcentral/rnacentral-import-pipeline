# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

import psycopg2
from pypika import Query, Table
from pypika import functions as fn


def balanced_ranges(counts, max_urs, max_accessions):
    ranges = []
    start = 1
    last = None
    urs = accessions = 0
    for urs_id, count in counts:
        if urs and (urs + 1 > max_urs or accessions + count > max_accessions):
            ranges.append((start, last))
            start = last + 1
            urs = accessions = 0
        urs += 1
        accessions += count
        last = urs_id

    if last is None:
        raise ValueError("No accessions to build ranges from")

    ranges.append((start, last))
    return ranges


def upi_taxid_ranges(ranges, tablename, db_url):
    table = Table(tablename)
    with psycopg2.connect(db_url) as conn:
        with conn.cursor() as cur:
            for start, stop in ranges:
                query = (
                    Query.from_(table)
                    .select(fn.Min(table.id), fn.Max(table.id))
                    .where(table.precompute_urs_id.between(start, stop))
                )

                cur.execute(str(query), ((start, stop),))
                results = cur.fetchall()
                assert len(results) == 1
                (upi_start, upi_max) = results[0]
                ## This can end up with nulls, catch it and skip writing
                ## here to avoid crashing the pipeline
                if upi_start is None:
                    continue
                yield (start, upi_start, upi_max)


def write(handle, output, tablename="precompute_urs_taxid", db_url=None):
    reader = csv.reader(handle)
    ranges = [(start, stop) for (_tablename, start, stop) in reader]
    writer = csv.writer(output)
    writer.writerows(upi_taxid_ranges(ranges, tablename, db_url))
