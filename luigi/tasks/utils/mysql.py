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

import urlparse

from plumbum import local

from tasks.utils.files import atomic_output


def download(uri, path):
    """
    Download the specific table using MySQL to a specific file. The input URI
    should be of the form:
        mysql://rfamro@mysql-rfam-public.ebi.ac.uk:4497/Rfam?database_link

    that is:
        mysql://username@host:port/database?table
    """

    parts = urlparse.urlparse(uri)
    database = parts.path[1:]

    mysql = local['mysql']
    query = 'select * from %s' % parts.query
    with atomic_output(path) as out:
        cmd = mysql[
            '--host', parts.hostname,
            '--port', str(parts.port),
            '--user', parts.username,
            '--database', database,
        ]

        ((cmd << query) > out)()
