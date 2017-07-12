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

from tasks.rfam.clans_csv import RfamClansCSV

CONTROL_FILE = """LOAD CSV
FROM '{filename}'
HAVING FIELDS
(
    rfam_clan_id,
    name,
    description,
    family_count
)
INTO {db_url}
TARGET COLUMNS
(
    rfam_clan_id,
    name,
    description,
    family_count
)
WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_rfam_clans (
    rfam_clan_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    name character varying(40) COLLATE pg_catalog."default" NOT NULL,
    description character varying(1000) COLLATE pg_catalog."default" NOT NULL,
    family_count integer NOT NULL
);
$$,
$$
truncate table load_rfam_clans;
$$

AFTER LOAD DO
$$ insert into rfam_clans (
    rfam_clan_id,
    name,
    description,
    family_count
) (
select
    rfam_clan_id,
    name,
    description,
    family_count
from load_rfam_clans
)
ON CONFLICT (rfam_clan_id) DO UPDATE SET
    name = excluded.name,
    description = excluded.description,
    family_count = excluded.family_count
$$,
$$
truncate table load_rfam_clans;
drop table load_rfam_clans;
$$
;
"""


class RfamPGLoadClans(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the Rfam clan CSV file. The import procedure will
    update any existing clans and will not produce duplicates.
    """
    def requires(self):
        return [
            RfamClansCSV(),
        ]

    def control_file(self):
        filename = RfamClansCSV().output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='load_rfam_clans')
        )
