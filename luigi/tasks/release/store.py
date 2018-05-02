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
from rnacentral.db import get_db_connection


class DatabaseUpdater(luigi.Task):  # pylint: disable=R0904
    """
    Base class for running pgsql procedures.
    """

    def output(self):
        """
        Check that a sentinel file exists.
        """
        filename = '%s.txt' % self.__class__.__name__
        return luigi.LocalTarget(os.path.join(output().base, 'cmds', filename))

    def create_sentinel_file(self):
        """
        Create a sentinel file on completion.
        To rerun the import, the sentinel file must be deleted.
        """
        with open(self.output().fn, 'w') as sentinel_file:
            sentinel_file.write('Done')


class UpdateAccessions(DatabaseUpdater):  # pylint: disable=R0904
    """
    Merge accessions data from loading table into main table.
    """

    @property
    def conn(self):
        if not hasattr(self, "_conn"):
            self._conn = get_db_connection(db())
        return self._conn

    def run(self):
        cur = self.conn.cursor()
        cur.execute("SET work_mem TO '256MB'")
        cur.execute('SELECT rnc_update.update_rnc_accessions()')
        self.create_sentinel_file()
        self.conn.close()


class UpdateReferences(DatabaseUpdater):  # pylint: disable=R0904
    """
    Merge literature references from loading table into main table.
    """

    @property
    def conn(self):
        if not hasattr(self, "_conn"):
            self._conn = get_db_connection(db())
        return self._conn

    def run(self):
        cur = self.conn.cursor()
        cur.execute("SET work_mem TO '256MB'")
        cur.execute('SELECT rnc_update.update_literature_references()')
        self.create_sentinel_file()
        self.conn.close()


class UpdateCoordinates(DatabaseUpdater):  # pylint: disable=R0904
    """
    Merge coordinates from loading table into main table.
    """

    @property
    def conn(self):
        if not hasattr(self, "_conn"):
            self._conn = get_db_connection(db())
        return self._conn

    sql = """
    INSERT INTO rnacen.rnc_coordinates AS t1 (
      accession,
      name,
      local_start,
      local_end,
      strand,
      assembly_id,
      id
    )
    SELECT
      load.accession,
      load.name,
      load.local_start,
      load.local_end,
      load.strand,
      assembly.assembly_id,
      NEXTVAL('rnc_coordinates_pk_seq')
    FROM rnacen.load_rnc_coordinates as load
    join ensembl_assembly assembly
    on
        assembly.assembly_id = load.assembly_id
    WHERE
      load.name is not null
    ON CONFLICT (accession, name, local_start, local_end, assembly_id)
    DO NOTHING
    """

    def run(self):
        cur = self.conn.cursor()
        cur.execute("SET work_mem TO '256MB'")
        cur.execute(self.sql)
        self.create_sentinel_file()
        self.conn.close()


class RunRelease(DatabaseUpdater):  # pylint: disable=R0904
    """
    Import data from the temporary load_* tables into the main tables.
    """

    @property
    def conn(self):
        if not hasattr(self, "_conn"):
            self._conn = get_db_connection(db())
        return self._conn

    sql = """
    SELECT dbid, id
    FROM rnacen.rnc_release
    WHERE status = 'L'
    ORDER BY id
    """

    create_index_sql = """
    CREATE INDEX IF NOT EXISTS load_rnacentral_all$database
    ON rnacen.load_rnacentral_all(database)
    """

    def requires(self):
        yield UpdateAccessions()
        yield UpdateReferences()
        yield UpdateCoordinates()

    def run(self):
        cur = self.conn.cursor()
        cur.execute("SET work_mem TO '256MB'")
        cur.execute(self.create_index_sql)
        cur.execute("SELECT rnc_update.prepare_releases('F')")
        cur.execute(self.sql)
        for result in cur.fetchall():
            print "INFO: Executing release %i from database %i" % (result[1], result[0])
            cur.execute('SELECT rnc_update.new_update_release(%s, %s)',
                        (result[0], result[1]))
            print "INFO: Committing..."
            self.conn.commit()
        self.create_sentinel_file()
        self.conn.close()
