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

from tasks.rfam.families_csv import RfamFamiliesCSV
from tasks.rfam.pgload_clans import RfamPGLoadClans

CONTROL_FILE = """LOAD CSV
FROM '{filename}' WITH ENCODING ISO-8859-14
HAVING FIELDS
(
    rfam_model_id,
    short_name,
    long_name,
    description [null if blanks],
    rfam_clan_id [null if blanks],
    seed_count,
    full_count,
    length,
    domain [null if blanks],
    is_suppressed,
    rna_type
)
INTO {db_url}
TARGET COLUMNS
(
    rfam_model_id,
    short_name,
    long_name,
    description,
    rfam_clan_id,
    seed_count,
    full_count,
    length,
    domain,
    is_suppressed,
    rna_type
)
SET
    search_path = '{search_path}'

WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_rfam_models (
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    long_name character varying(200) COLLATE pg_catalog."default" NOT NULL,
    description character varying(2000) COLLATE pg_catalog."default",
    seed_count integer NOT NULL,
    full_count integer NOT NULL,
    length integer NOT NULL,
    is_suppressed boolean NOT NULL,
    rfam_clan_id character varying(20) COLLATE pg_catalog."default",
    domain character varying(50) COLLATE pg_catalog."default",
    rna_type character varying(250) COLLATE pg_catalog."default",
    short_name character varying(50) COLLATE pg_catalog."default"
);
$$,
$$
truncate table load_rfam_clans;
$$

AFTER LOAD DO
$$ insert into rfam_models (
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type
) (
select
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type
from load_rfam_models
)
ON CONFLICT (rfam_model_id) DO UPDATE SET
    short_name = excluded.short_name,
    long_name = excluded.long_name,
    description = excluded.description,
    seed_count = excluded.seed_count,
    full_count = excluded.full_count,
    length = excluded.length,
    is_suppressed = excluded.is_suppressed,
    rfam_clan_id = excluded.rfam_clan_id,
    domain = excluded.domain,
    rna_type = excluded.rna_type
$$,
$$
truncate table load_rfam_models;
$$,
$$
drop table load_rfam_models;
$$
;
"""


class RfamPGLoadFamilies(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the Rfam family CSV file. The importing will
    update any existing families and will not produce duplicates.
    """

    def requires(self):
        return [
            RfamFamiliesCSV(),
            RfamPGLoadClans(),
        ]

    def control_file(self):
        filename = RfamFamiliesCSV().output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='load_rfam_models'),
            search_path=self.db_search_path(),
        )
