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

import os

import luigi

from tasks.config import output
from tasks.utils.pgloader import PGLoader

from .generic import file_pattern


CONTROL_FILE = """
LOAD CSV
FROM ALL FILENAMES MATCHING ~<{pattern}>
IN DIRECTORY '{directory}'
HAVING FIELDS (
    CRC64,
    LEN,
    {sequence_column},
    DATABASE,
    AC,
    OPTIONAL_ID,
    VERSION,
    TAXID,
    MD5
)
INTO {db_url}
TARGET COLUMNS (
    CRC64,
    LEN,
    {sequence_column},
    DATABASE,
    AC,
    OPTIONAL_ID,
    VERSION,
    TAXID,
    MD5
)

WITH
    batch rows = 25000,
    batch size =  512MB,
    workers = 10,
    concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    search_path = '{search_path}',
    work_mem to '256 MB',
    maintenance_work_mem to '1 GB'
;
"""


class PGLoadSequences(PGLoader):  # pylint: disable=R0921,R0904
    """
    This is the base class for loading sequences (short or long) into the
    database. This will only populate the load_rnacentral_all table, and will
    not run further scripts.
    """

    database = luigi.Parameter(default='all')
    type = luigi.ChoiceParameter(choice=['short', 'long'])

    def sequence_column(self):
        """
        This is the name of the sequence column. This is computed from the name
        of the directory, which should be one of short or long.
        """
        return 'SEQ_' + self.type.upper()

    def control_file(self):
        return CONTROL_FILE.format(
            pattern=file_pattern(self.database),
            db_url=self.db_url(table='load_rnacentral_all'),
            search_path=self.db_search_path(),
            directory=os.path.join(output().base, self.type),
            sequence_column=self.sequence_column(),
        )
