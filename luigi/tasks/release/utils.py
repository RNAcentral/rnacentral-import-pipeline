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

from contextlib import contextmanager

import psycopg2

from tasks.config import db


class DatabaseConnection(object):
    """
    Some general utilities for having a connection to the database.
    """

    @contextmanager
    def connection(self):
        """
        Opens a connection to the datbase.
        """

        connection = psycopg2.connect(db().psycopg2_string())
        connection.set_session(autocommit=False)
        yield connection
        connection.close()

    @contextmanager
    def cursor(self, commit_on_leave=True):
        """
        Create a cursor to the database. The cursor will be built with the
        standard options for work_mem and maintenance_work_mem. If
        commit_on_leave is true then the cursor will commit once the context
        handler exits.
        """

        with self.connection() as connection:
            cursor = connection.cursor()
            cursor.execute("SET work_mem to '256MB'")
            cursor.execute("SET maintenance_work_mem TO '512MB'")
            yield cursor
            if commit_on_leave:
                cursor.commit()
            cursor.close()
