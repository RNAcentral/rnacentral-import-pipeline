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

from tasks.ontologies import Ontologies

from .quickgo_data import QuickGoData


CONTROL_FILE = """
LOAD CSV
FROM '{filename}'
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
    rna_id varchar(50),
    qualifier varchar(30),
    assigned_by varchar(50),
    extensions jsonb,
    ontology_term_id varchar(15),
    evidence_code varchar(15)
);
$$,
$$
truncate table {load_table};
$$

AFTER LOAD DO
$$
INSERT INTO {final_table} (
    rna_id,
    qualifier,
    assigned_by,
    extensions,
    ontology_term_id,
    evidence_code
) (
SELECT
    rna_id,
    qualifier,
    assigned_by,
    extensions,
    ontology_term_id,
    evidence_code
FROM {load_table}
)
ON CONFLICT (rna_id, qualifier, assigned_by, ontology_term_id, evidence_code)
DO UPDATE
SET
    extensions = excluded.extensions
;
$$,
$$
DROP TABLE {load_table};
$$
;
"""


class QuickGoLoadAnnotations(PGLoader):
    def requires(self):
        return [
            QuickGoData(),
            Ontologies(),
        ]

    def control_file(self):
        output = self.requires()[0].output()
        table = 'go_term_annotations'
        load_table = 'load_' + table
        fields = ', '.join(output.annotations.headers)

        return CONTROL_FILE.format(
            filename=output.annotations.filename,
            fields=fields,
            columns=fields,
            final_table=table,
            load_table=load_table,
            db_url=self.db_url(table=load_table),
            search_path=self.db_search_path(),
        )
