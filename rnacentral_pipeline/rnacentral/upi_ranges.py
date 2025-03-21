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

from ..db import cursor


def ranges_between(start, stop, max_size):
    last = None
    for start in range(start, stop, max_size):
        last = min(start + max_size, stop)
        yield (start, last)

    if last is None:
        raise ValueError("Some iterated over nothing")

    if last < stop:
        yield (last, stop)


def upi_ranges(dbconf, table_name, max_size):
    """
    This will create range of the ids for all UPI's in the database.
    """

    with cursor(dbconf) as cur:
        cur.execute("select max(id) from %s" % table_name)
        stop = cur.fetchone()[0]

    return ranges_between(1, stop, max_size)


def to_file(dbconf, table_name, max_size, output):
    ranges = upi_ranges(dbconf, table_name, max_size)
    ranges = [[table_name] + list(r) for r in ranges]
    csv.writer(output).writerows(ranges)
