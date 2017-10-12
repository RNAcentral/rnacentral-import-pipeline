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

import os

import luigi

from tasks.config import db, output
from .utils.db import get_db_connection


SQL = """
select
    dbid,
    id
from rnacen.rnc_release
where status = 'L'
order by id
"""

CREATE_INDEX = """
create index if not exists load_rnacentral_all$database
on rnacen.load_rnacentral_all(database)
"""


class StoreRelease(luigi.Task):  # pylint: disable=R0904
    """
    Import data from the temporary load_* tables into the main tables.
    Create a sentinel file on completion.
    To rerun the import, the sentinel file must be deleted.
    """

    def run(self):
        conn = get_db_connection(db())
        cur = conn.cursor()
        cur.execute("SET work_mem TO '256MB'")
        cur.execute(CREATE_INDEX)
        cur.execute("SELECT rnc_update.prepare_releases('F')")
        cur.execute(SQL)
        for result in cur.fetchall():
            print "Running release %i from database %s" % (result[1], result[0])
            cur.execute('SELECT rnc_update.new_update_release(%s, %s)',
                        (result[0], result[1]))
            print "Committing..."
            conn.commit()
        cur.execute('SELECT rnc_update.update_rnc_accessions()')
        cur.execute('SELECT rnc_update.update_literature_references()')
        conn.close()
        with open(self.output().fn , 'w') as sentinel_file:
            sentinel_file.write('Done')


    def output(self):
        """
        Check that a sentinel file exists.
        """
        return luigi.LocalTarget(os.path.join(output().base, '%s.txt' % self.__class__.__name__))
