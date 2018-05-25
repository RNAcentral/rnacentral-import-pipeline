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

import psycopg2
import psycopg2.extras

from tasks.config import db


def count(sql):
    """
    Returns the count of the nubmer of elements returned by the given sql
    statement. This will open a new connection to the database.
    """

    connection = psycopg2.connect(db().pgloader_url())
    cursor = connection.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute('select count(*) total from (%s) t' % sql)
    result = cursor.fetchone()
    return result['total']
