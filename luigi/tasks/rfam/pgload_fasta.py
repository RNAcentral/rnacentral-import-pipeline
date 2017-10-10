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

from .fasta_csv import RfamFastaCSV

CONTROL_FILE = """LOAD CSV
FROM '{filename}'
HAVING FIELDS
(
    upi,
    date
)
INTO {db_url}
TARGET COLUMNS
(
    upi,
    date
)
WITH
    skip header = 1,
    fields escaped by double-quote,
    fields terminated by ','
SET
    search_path = '{search_path}'
BEFORE LOAD DO
$$
create table if not exists load_rfam_analyzed_sequences
    upi character varying(13) COLLATE pg_catalog."default" NOT NULL,
    date date NOT NULL
);
$$,
$$
truncate table load_rfam_model_hits;
$$

AFTER LOAD DO
$$
insert into rfam_analyzed_sequences (
    upi,
    date,
) (
select
    upi,
    date
from load_rfam_analyzed_sequences
)
ON CONFLICT (upi) DO UPDATE SET
    upi = excluded.upi,
    date = excluded.date
$$,
$$
drop table load_rfam_analyzed_sequences;
$$
;
"""


class RfamPGLoadFasta(PGLoader):  # pylint: disable=R0904
    """
    This will run pgloader on the CSV file produced by FastaCSV. This fills out
    the rfam_analyzed_sequences table and the procedure will fail if the same
    sequence has been run more than once.
    """

    def requires(self):
        return [
            RfamFastaCSV(),
        ]

    def control_file(self):
        filename = RfamFastaCSV().output().fn
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table='load_rfam_analyzed_sequences'),
            search_path=self.db_search_path(),
        )
