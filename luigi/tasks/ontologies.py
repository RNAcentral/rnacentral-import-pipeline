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

from tasks.config import ontologies

from tasks.utils.pgloader import PGLoader

from tasks.rfam.go_terms import RfamGoTerms
from tasks.quickgo.quickgo_data import QuickGoData


CONTROL_FILE = """
LOAD CSV
FROM ALL FILENAMES MATCHING ~<{pattern}> WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    ontology_term_id,
    ontology,
    name,
    definition
)
INTO {db_url}
TARGET COLUMNS
(
    ontology_term_id,
    ontology,
    name,
    definition
)
SET
    search_path = '{search_path}'

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
evidence
create table if not exists {load_table} (
    ontology_term_id varchar(15),
    ontology varchar(5),
    name text,
    definition text
);
$$,
$$
truncate table {load_table};
$$

AFTER LOAD DO
$$ insert into {final_table} (
    ontology_term_id,
    ontology,
    name,
    definition
) (
select
    ontology_term_id,
    ontology,
    name,
    definition
from load_go_terms
)
ON CONFLICT (ontology_term_id) DO UPDATE SET
    ontology_term_id = excluded.ontology_term_id,
    ontology = excluded.ontology,
    name = excluded.name,
    definition = excluded.definition
;
$$,
$$
drop table {load_table};
$$
;
"""


class Ontologies(PGLoader):  # pylint: disable=R0904

    def requires(self):
        return [
            QuickGoData(),
            RfamGoTerms(),
        ]

    def control_file(self):
        table = 'ontology_terms'
        load_table = 'load_' + table
        return CONTROL_FILE.format(
            pattern=ontologies.to_load('*'),
            final_table=table,
            load_table=load_table,
            search_path=self.db_search_path(),
            db_url=self.db_url(table=load_table),
        )
