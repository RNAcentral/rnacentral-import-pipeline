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


class PGLoadAccessions(PGLoader):  # pylint: disable=R0921,R0904
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
        'ensembl.*.csv' will load all ensembl accessions.
        """
        pass

    def control_file(self):
        config = output()
        directory = os.path.join(config.base, 'ac_info')
        return CONTROL_FILE.format(
            pattern=self.pattern(),
            db_url=self.db_url(table='load_rnc_accessions'),
            search_path=self.db_search_path(),
            directory=directory,
        )
