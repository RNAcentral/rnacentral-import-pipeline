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

from tasks.publications.pubmed import PubmedLoader

from .quickgo_data import QuickGoData
from .load_annotations import QuickGoLoadAnnotations

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
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    pubmed_id int
);
$$,
$$
truncate table {load_table};
$$

AFTER LOAD DO
$$ insert into {final_table} (go_term_annotation_id, ref_pubmed_id)
(
select
    annotations.go_term_annotation_id,
    {load_table}.pubmed_id
from {load_table}
join go_term_annotations annotations
on
    annotations.rna_id = {load_table}.rna_id
    AND annotations.qualifier = {load_table}.qualifier
    AND annotations.assigned_by = {load_table}.assigned_by
    AND annotations.ontology_term_id = {load_table}.ontology_term_id
    AND annotations.evidence_code = {load_table}.evidence_code
)
ON CONFLICT (go_term_annotation_id, ref_pubmed_id)
DO NOTHING
;
$$,
$$
drop table {load_table};
$$
;
"""


class QuickGoLoadPublicationMapping(PGLoader):
    def requires(self):
        return [
            QuickGoData(),
            PubmedLoader(),
            QuickGoLoadAnnotations(),
        ]

    def control_file(self):
        output = self.requires()[0].output()
        table = 'go_term_publication_map'
        load_table = 'load_' + table
        fields = ', '.join(output.publication_mappings.headers)

        return CONTROL_FILE.format(
            filename=output.publication_mappings.filename,
            directory=publications().to_load(),
            final_table=table,
            load_table=load_table,
            db_url=self.db_url(table=load_table),
            columns=fields,
            fields=fields,
            search_path=self.db_search_path(),
        )
