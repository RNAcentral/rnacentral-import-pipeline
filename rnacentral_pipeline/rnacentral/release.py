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
