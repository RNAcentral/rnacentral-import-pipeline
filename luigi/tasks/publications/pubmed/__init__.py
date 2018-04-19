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

from tasks.config import publications

from tasks.utils.pgloader import PGLoader

from tasks.quickgo.quickgo_data import QuickGoData


CONTROL_FILE = """
LOAD CSV
FROM ALL FILENAMES MATCHING ~<{pattern}>
IN DIRECTORY '{directory}'
WITH ENCODING ISO-8859-14

HAVING FIELDS ({fields})
INTO {db_url}
TARGET COLUMNS ({columns})

SET
    search_path = '{search_path}'

WITH
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists {load_table} (
    ref_pubmed_id int,
    authors text,
    location varchar(4000),
    title varchar(4000)
);
$$,
$$
truncate table {load_table};
$$

AFTER LOAD DO
$$ insert into {final_table} (
    ref_pubmed_id,
    authors,
    location,
    title
) (
select distinct
    ref_pubmed_id,
    authors,
    location,
    title
from {load_table}
)
ON CONFLICT (ref_pubmed_id) DO UPDATE SET
    ref_pubmed_id = excluded.ref_pubmed_id,
    authors = excluded.authors,
    location = excluded.location,
    title = excluded.title
;
$$,
$$
drop table {load_table};
$$
;
"""


class PubmedLoader(PGLoader):
    def requires(self):
        return [
            QuickGoData(),
        ]

    def control_file(self):
        output = self.requires()[0].output()
        table = 'ref_pubmed'
        load_table = 'load_' + table
        fields = ','.join(output.publications.header)

        return CONTROL_FILE.format(
            pattern='.*csv',
            directory=publications().to_load(),
            final_table=table,
            load_table=load_table,
            db_url=self.db_url(table=load_table),
            columns=fields,
            fields=fields,
            search_path=self.db_search_path(),
        )
