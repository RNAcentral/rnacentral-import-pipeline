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
from pypika import Field, Query, Table
from pypika import functions as fn


def upi_taxid_ranges(ranges, tablename, db_url):
    table = Table(tablename)
    with psycopg2.connect(db_url) as conn:
        with conn.cursor() as cur:
            for (start, stop) in ranges:
                query = (
                    Query.from_(table)
                    .select(fn.Min(table.id), fn.Max(table.id))
                    .where(table.precompute_urs_id.between(start, stop))
                )

                cur.execute(str(query), ((start, stop),))
                results = cur.fetchall()
                assert len(results) == 1
                (upi_start, upi_max) = results[0]
                yield (start, upi_start, upi_max)


def write(handle, output, tablename="precompute_urs_taxid", db_url=None):
    reader = csv.reader(handle)
    ranges = [(start, stop) for (_tablename, start, stop) in reader]
    writer = csv.writer(output)
    writer.writerows(upi_taxid_ranges(ranges, tablename, db_url))
