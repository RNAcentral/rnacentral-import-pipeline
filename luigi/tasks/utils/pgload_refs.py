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
import abc

from tasks.utils.pgloader import PGLoader
from tasks.config import output


CONTROL_FILE = """
LOAD CSV
FROM all filenames matching ~<{pattern}>
in directory '{directory}'
HAVING FIELDS  (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi
)
INTO {db_url}
TARGET COLUMNS (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi
)

WITH
    truncate,
    batch rows = 500,
    batch size = 32MB,
    prefetch rows = 500,
    workers = 4, concurrency = 2,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

SET
    work_mem to '256 MB',
    maintenance_work_mem to '256 MB',
    search_path = '{search_path}'

BEFORE LOAD DO
$$
ALTER TABLE rnacen.load_rnc_references SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_references SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$
;
"""


class PGLoadReferences(PGLoader):  # pylint: disable=R0921,R0904
    """
    This is a base class to load accessions. By specifying the pattern this can
    be used to load all accessions, or only those for a specific database.
    """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def pattern(self):
        """
        This should give the pattern of filenames to load. If this returns
        '.*.csv' then it iwll load all accessions. Using something like
        'ensembl.*.csv' will load all ensembl references.
        """
        pass

    def control_file(self):
        config = output()
        directory = os.path.join(config.base, 'refs')
        return CONTROL_FILE.format(
            pattern=self.pattern(),
            db_url=self.db_url(table='load_rnc_references'),
            search_path=self.db_search_path(),
            directory=directory,
        )
