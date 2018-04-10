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


from tasks.config import db
from .store import UpdateCoordinates
from .store import DatabaseUpdater
from rnacentral.db import get_db_connection

SQL = """
update rnc_coordinates coord
set
    "name" = coord.primary_accession,
    primary_start = coord.local_start,
    primary_end = coord.local_end
from rnc_accessions acc
where
    acc.accession = coord.accession
    and acc."database" in ('ENSEMBL', 'GENCODE')
"""


class EnsemblCoordinates(DatabaseUpdater):  # pylint: disable=R0904
    """
    This will update all Ensembl coordinates in the database. Basically we can
    copy part of the Ensembl/GENCODE data to get complete information about the
    coordinates, no need for external API's.
    """

    def requires(self):
        yield UpdateCoordinates()

    @property
    def conn(self):
        """
        Get a connection to the database.
        """

        if not hasattr(self, "_conn"):
            self._conn = get_db_connection(db())  # pylint: disable=W0201
        return self._conn

    def run(self):
        cur = self.conn.cursor()
        cur.execute(SQL)
        self.conn.close()
        self.create_sentinel_file()
