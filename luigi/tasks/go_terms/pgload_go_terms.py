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

from tasks.utils.pgloader import PGLoader

from tasks.go_terms.go_terms_csv import GoTermsCSV

CONTROL_FILE = """LOAD CSV
FROM '{filename}' WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    go_term_id,
    name
)
INTO {db_url}
TARGET COLUMNS
(
    go_term_id,
    name
)
SET
    search_path = '{search_path}'

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_go_terms (
    go_term_id varchar(10) NOT NULL,
    name text COLLATE pg_catalog."default" NOT NULL
);
$$,
$$
truncate table load_go_terms;
$$

AFTER LOAD DO
$$ insert into go_terms (
    go_term_id,
    name
) (
select
    go_term_id,
    name
from load_go_terms
)
ON CONFLICT (go_term_id) DO UPDATE SET
    go_term_id = excluded.go_term_id,
    name = excluded.name
;
$$,
$$
drop table load_go_terms;
$$
;
"""


class PGLoadGoTerms(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the Rfam family CSV file. The importing will
    update any existing families and will not produce duplicates.
    """

    def requires(self):
        return [
            GoTermsCSV(),
        ]

    def control_file(self):
        filename = GoTermsCSV().output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='load_go_terms'),
            search_path=self.db_search_path(),
        )
