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

from tasks.go_terms.pgload_go_terms import PGLoadGoTerms
from .go_term_mapping_csv import RfamGoTermsCSV
from .pgload_families import RfamPGLoadFamilies

CONTROL_FILE = """LOAD CSV
FROM '{filename}' WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    go_term_id,
    rfam_model_id
)
INTO {db_url}
TARGET COLUMNS
(
    go_term_id,
    rfam_model_id
)
SET
    search_path = '{search_path}'

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_rfam_go_terms (
    go_term_id character varying(10) COLLATE pg_catalog."default" NOT NULL,
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL
);
$$,
$$
truncate table load_rfam_go_terms;
$$

AFTER LOAD DO
$$ insert into rfam_go_terms (
    go_term_id,
    rfam_model_id
) (
select
    go_term_id,
    rfam_model_id
from load_rfam_go_terms
)
ON CONFLICT (go_term_id, rfam_model_id) DO UPDATE SET
    go_term_id = excluded.go_term_id,
    rfam_model_id = excluded.rfam_model_id
;
$$,
$$
drop table load_rfam_go_terms;
$$
;
"""


class RfamPGLoadGoTerms(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the Rfam go term mapping  CSV file. The importing
    will update any existing mappings and will not produce duplicates.
    """

    def requires(self):
        return [
            RfamGoTermsCSV(),
            PGLoadGoTerms(),
            RfamPGLoadFamilies(),
        ]

    def control_file(self):
        filename = RfamGoTermsCSV().output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='load_rfam_go_terms'),
            search_path=self.db_search_path(),
        )
