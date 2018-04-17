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

from .eco_code_csv import EcoCodeCSV


CONTROL_FILE = """LOAD CSV
FROM '{filename}' WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    eco_term_id,
    name,
    description
)
INTO {db_url}
TARGET COLUMNS
(
    eco_term_id,
    name,
    description
)
SET
    search_path = '{search_path}'

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists {table} (
    eco_term_id varchar(10) NOT NULL,
    name text COLLATE pg_catalog."default" NOT NULL,
    description text
);
$$,
$$
truncate table {table};
$$

AFTER LOAD DO
$$ insert into {table} (
    eco_term_id,
    name,
    description
) (
select
    eco_term_id,
    name,
    description
from {table}
)
ON CONFLICT (eco_term_id) DO UPDATE SET
    eco_term_id = excluded.eco_term_id,
    name = excluded.name,
    description = excluded.description
;
$$,
$$
drop table {table};
$$
;
"""


class PGLoadEcoTerms(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the Rfam family CSV file. The importing will
    update any existing families and will not produce duplicates.
    """

    def requires(self):
        return EcoCodeCSV()

    def control_file(self):
        filename = EcoCodeCSV().output().fn
        table = 'load_eco_codes'
        return CONTROL_FILE.format(
            filename=filename,
            table=table,
            db_url=self.db_url(table=table),
            search_path=self.db_search_path(),
        )
