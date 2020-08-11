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

import json
import logging

from ..db import connection

LOGGER = logging.getLogger(__name__)


CREATE_INDEX_SQL = """
CREATE INDEX IF NOT EXISTS load_rnacentral_all$database
ON rnacen.load_rnacentral_all(database)
"""

TO_RELEASE = """
SELECT dbid, id
FROM rnacen.rnc_release
WHERE status = 'L'
ORDER BY id
"""

COUNT_QUERY = """
SELECT
    db.descr,
    count(distinct xref.upi)
from xref
join rnc_database db 
on
    db.id = xref.dbid
where
    xref.deleted = 'N'
group by db.descr
"""

LOAD_COUNT_QUERY = """
SELECT
    load.database,
    count(distinct load.md5)
from load_rnacentral load
group by database
"""


def run(db_url):
    """
    Run the release logic. Basically this will run the select commands that are
    needed to update all the tables and then run the release update logic.
    """

    with connection(db_url) as conn:
        cursor = conn.cursor()
        cursor.execute("SET work_mem TO '256MB'")
        cursor.execute('SELECT rnc_update.update_rnc_accessions()')
        cursor.execute('SELECT rnc_update.update_literature_references()')
        cursor.execute(CREATE_INDEX_SQL)
        cursor.execute("SELECT rnc_update.prepare_releases('F')")
        cursor.execute(TO_RELEASE)
        for (dbid, rid) in cursor.fetchall():
            LOGGER.info("Executing release %i from database %i", rid, dbid)
            cursor.execute('SELECT rnc_update.new_update_release(%s, %s)',
                           (dbid, rid))
            conn.commit()


def check(limit_file, db_url, default_allowed_change=0.30):
    """
    Check the load tables for reasonable looking sequence counts.
    """

    limits = json.load(limit_file)
    cur_counts = {}
    new_counts = {}
    with connection(db_url) as conn:
        cursor = conn.cursor()
        cursor.execute(COUNT_QUERY)
        for (descr, raw_count) in cursor.fetchall():
            cur_counts[descr] = float(raw_count)

        cursor.execute(LOAD_COUNT_QUERY)
        for (descr, raw_count) in cursor.fetchall():
            new_counts[descr] = float(raw_count)

    problems = False
    for name, previous in cur_counts.items():
        current = new_counts[name]
        change = (current - previous) / float(current)
        if change > limits.get(name, default_allowed_change):
            LOGGER.error("Database %s increased by %f", name, change)
            problems = True

    if problems:
        raise ValueError("Overly large changes with release")
