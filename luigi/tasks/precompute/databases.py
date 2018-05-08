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

from hashlib import sha256
from luigi.contrib.postgres import PostgresQuery

import luigi

from tests.tasks.helpers import count
from tasks.config import db


SQL = """
WITH t1 AS (
	SELECT
	xref.upi,
	xref.taxid,
	STRING_AGG(DISTINCT display_name, ',' ORDER BY display_name) AS databases
	FROM xref, rnc_database, rna
	WHERE
		deleted = 'N'
		AND xref.dbid=rnc_database.id
		AND rna.upi = xref.upi
		AND rna.id > {min_id} AND rna.id <= {max_id}
	GROUP BY xref.upi, taxid
)
UPDATE rnc_rna_precomputed t2
SET databases = t1.databases
FROM t1
WHERE
	t1.upi = t2.upi
	AND t2.taxid = t2.taxid
"""


class PrecomputeDatabaseList(PostgresQuery):
	query = luigi.Parameter()
	update_id = luigi.Parameter()
	host = db().host
	database = db().db_name
	user = db().user
	password = db().password
	table = 'rnc_rna_precomputed__databases'


class PrecomputeDatabaseListWrapper(luigi.WrapperTask):
	min_id = luigi.IntParameter(default=0)
	max_id = luigi.IntParameter(default=0)
	chunk_size = luigi.IntParameter(default=1000)

	def requires(self):
		i = self.min_id
		if self.max_id == 0:
			self.max_id = count('SELECT * FROM rna')
		while i < self.max_id:
			query = SQL.format(min_id=i, max_id=i+self.chunk_size)
			yield PrecomputeDatabaseList(query=query,
			                             update_id=sha256(query).hexdigest())
			i += self.chunk_size
