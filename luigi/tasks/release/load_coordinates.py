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

from tasks.config import output
from tasks.utils.pgloader import PGLoader
from .manage_files import SplitMergedFile


CONTROL_FILE = """
LOAD CSV
FROM ALL FILENAMES MATCHING ~<{pattern}>
IN DIRECTORY '{directory}'
HAVING FIELDS (
    accession,
    primary_accession,
    local_start,
    local_end,
    strand
)
INTO {db_url}
TARGET COLUMNS (
    accession,
    primary_accession,
    local_start,
    local_end,
    strand
)

WITH truncate,
    drop indexes,
    batch rows = 500,
    batch size = 32MB,
    prefetch rows = 500,
    workers = 2, concurrency = 1,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 GB',
    search_path = '{search_path}'

BEFORE LOAD DO
$$
ALTER TABLE rnacen.load_rnc_coordinates SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_coordinates SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$
;
"""


class LoadCoordinates(PGLoader):  # pylint: disable=R0904
    """
    This will load coordinates. The database parameter defaults to all
    coordinates, if a value is given then it is assumed to be the name of the
    database to load. All files that begin with that name will be loaded.
    """
    directory = 'genomic_locations'

    def requires(self):
        return SplitMergedFile(directory=self.directory)

    def control_file(self):
        config = output()
        directory = os.path.join(config.base, self.directory)
        return CONTROL_FILE.format(
            pattern='chunk_*',
            db_url=self.db_url(table='load_rnc_coordinates'),
            search_path=self.db_search_path(),
            directory=directory,
        )
