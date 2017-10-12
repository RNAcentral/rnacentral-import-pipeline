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
from .utils.generic import file_pattern
from .manage_files import SplitFiles


CONTROL_FILE = """
LOAD CSV
FROM ALL FILENAMES MATCHING ~<{pattern}>
IN DIRECTORY '{directory}'
HAVING FIELDS (
    accession,
    parent_ac,
    seq_version,
    feature_start,
    feature_end,
    feature_name,
    ordinal [null if ""],
    is_composite,
    non_coding_id,
    database,
    external_id,
    optional_id,
    project,
    division,
    keywords,
    description,
    species,
    common_name,
    organelle,
    classification,
    allele,
    anticodon,
    chromosome,
    experiment,
    function,
    gene,
    gene_synonym,
    inference,
    locus_tag,
    map,
    mol_type,
    ncrna_class,
    note,
    old_locus_tag,
    operon,
    product,
    pseudogene,
    standard_name,
    db_xref
)
INTO {db_url}
TARGET COLUMNS (
    accession,
    parent_ac,
    seq_version,
    feature_start,
    feature_end,
    feature_name,
    ordinal,
    is_composite,
    non_coding_id,
    database,
    external_id,
    optional_id,
    project,
    division,
    keywords,
    description,
    species,
    common_name,
    organelle,
    classification,
    allele,
    anticodon,
    chromosome,
    experiment,
    function,
    gene,
    gene_synonym,
    inference,
    locus_tag,
    map,
    mol_type,
    ncrna_class,
    note,
    old_locus_tag,
    operon,
    product,
    pseudogene,
    standard_name,
    db_xref
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
ALTER TABLE rnacen.load_rnc_accessions SET (
    autovacuum_enabled = false,
    toast.autovacuum_enabled = false
);
$$

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_rnc_accessions SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$
;
"""


class LoadAccessions(PGLoader):  # pylint: disable=R0904
    """
    This will load accessions. The database parameter defaults to all
    acecssion, if a value is given then it is assumed to be the name of the
    database to load. All files that begin with that name will be loaded.
    """
    database = luigi.Parameter(default='all')
    directory = 'ac_info'

    def requires(self):
        return SplitFiles(directory=self.directory)

    def control_file(self):
        config = output()
        directory = os.path.join(config.base, self.directory)
        return CONTROL_FILE.format(
            pattern=file_pattern(self.database),
            db_url=self.db_url(table='load_rnc_accessions'),
            search_path=self.db_search_path(),
            directory=directory,
        )
