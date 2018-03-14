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
from luigi.local_target import atomic_file

from rnacentral.db import get_db_connection
from tasks.config import db
from tasks.config import export

CREATE_ACTIVE_TABLE = """
CREATE TEMP TABLE rna_active (
    upi varchar(30) primary KEY
)
"""

INSERT_INTO_ACTIVE_TABLE = """
INSERT INTO rna_active
select distinct upi from xref_p{dbid}_not_deleted
ON CONFLICT DO NOTHING
"""

QUERY = """
select
    hits.upi,
    hits.rfam_model_id,
    score,
    e_value,
    sequence_start,
    sequence_stop,
    model_start,
    model_stop,
    models.long_name
from rfam_model_hits hits
join rna_active active on active.upi = hits.upi
join rfam_models models on models.rfam_model_id = hits.rfam_model_id
"""



class RfamAnnotations(luigi.Task):
    def output(self):
        filename = export().rfam('rfam_annotations.tsv')
        return luigi.LocalTarget(filename)

    def command(self, example=False):
        query = QUERY.replace('\n', '')
        if example:
            query = '%s LIMIT 10' % query
        return "COPY ({query}) to STDOUT DELIMITER '\t'".format(
            query=QUERY.replace('\n', ' '),
        )

    def populate_active_table(self, connection):
        cursor = connection.cursor()
        cursor.execute(CREATE_ACTIVE_TABLE)
        cursor.execute('select id from rnc_database order by id')
        dbids = [r[0] for r in cursor.fetchall()]
        for dbid in dbids:
            cursor.execute(INSERT_INTO_ACTIVE_TABLE.format(dbid=dbid))

    def run(self):
        connection = get_db_connection(db(), connect_timeout=10 * 60)
        # self.populate_active_table(connection)
        with atomic_file(self.output().fn) as out:
            cursor = connection.cursor()
            cursor.copy_expert(self.command(), out)
            connection.close()
