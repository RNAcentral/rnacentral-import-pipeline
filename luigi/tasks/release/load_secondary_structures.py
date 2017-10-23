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
import logging
from glob import glob

import luigi

from tasks.config import output
from tasks.utils.pgloader import PGLoader

from .manage_files import SplitMergedFile
from .utils.generic import file_pattern


LOGGER = logging.getLogger(__name__)


CONTROL_FILE = """
LOAD CSV
FROM
    ALL FILENAMES
    MATCHING ~<{pattern}>
    IN DIRECTORY '{directory}'
HAVING FIELDS (
    rnc_accession_id,
    secondary_structure,
    md5
) INTO {db_url}
TARGET COLUMNS (
    rnc_accession_id,
    secondary_structure,
    md5
)
WITH
    drop indexes,
    SKIP HEADER = 1,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','
SET
    search_path = '{search_path}'

BEFORE LOAD DO
$$
create table if not exists load_{tablename} (
    rnc_accession_id varchar(100),
    secondary_structure text,
    md5 varchar(32)
);
$$,
$$
truncate table load_{tablename};
$$

AFTER LOAD DO
$$
INSERT INTO {tablename} (
    rnc_accession_id,
    secondary_structure,
    md5
) (
select distinct
    rnc_accession_id,
    secondary_structure,
    md5
from load_{tablename}
)
-- We can't include two ON CONFLICT statements so this will do extra inserts
ON CONFLICT (rnc_accession_id, md5) DO UPDATE SET
    rnc_accession_id = excluded.rnc_accession_id,
    secondary_structure = excluded.secondary_structure,
    md5 = excluded.md5
;
$$,
$$
drop table load_{tablename};
$$
;
"""


class LoadSecondaryStructures(PGLoader):  # pylint: disable=R0904
    """
    This will load only the given file in the secondary structure directory.
    """

    database = luigi.Parameter(default='all')
    directory = 'secondary_structure'

    def requires(self):
        return SplitMergedFile(directory=self.directory)

    def data_directory(self):
        """
        Gets the path to the directory of secondary structures.
        """
        return os.path.join(output().base, self.directory)

    def file_pattern(self):
        """
        Build the pattern of filenames to use.
        """
        return file_pattern(self.database)

    def data_files(self):
        """
        Get the fill path, with pattern to the filenames to use.
        """
        return os.path.join(self.data_directory(), self.file_pattern())

    def control_file(self):
        tablename = 'rnc_secondary_structure'
        return CONTROL_FILE.format(
            pattern=self.file_pattern(),
            directory=self.data_directory(),
            tablename=tablename,
            db_url=self.db_url('load_%s' % tablename),
            search_path=self.db_search_path(),
        )

    def run(self):
        """
        It is possible that the directory of secondary structures does not
        exist. If this is the case then nothing is actually run.
        """

        if glob(self.data_files()):
            return super(LoadSecondaryStructures, self).run()
        LOGGER.info("No secondary structures found, skipping")
