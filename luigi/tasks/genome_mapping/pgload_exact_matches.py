# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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

from tasks.utils.pgloader import PGLoader
from .genome_mapping_tasks import ParsePslOutput
from tasks.config import genome_mapping


CONTROL_FILE = """LOAD CSV
FROM '{filename}'
HAVING FIELDS
(
    rna_id,
    upi,
    taxid,
    unused_field,
    region_id,
    chromosome,
    start,
    stop,
    strand,
    assembly_id,
    identity
)
INTO {db_url}
TARGET COLUMNS
(
    rna_id,
    upi,
    taxid,
    region_id,
    chromosome,
    start,
    stop,
    strand,
    assembly_id,
    identity
)
WITH
    fields terminated by '\t'
SET
    search_path = 'rnacen'
BEFORE LOAD DO
$$ DROP TABLE IF EXISTS {table}; $$,
$$ CREATE TABLE IF NOT EXISTS rnacen.{table} (
    rna_id varchar(50) NULL,
    upi varchar(13) NULL,
    taxid int4 NULL,
    region_id varchar(100) NULL,
    chromosome varchar(100) NULL,
    "start" int4 NULL,
    stop int4 NULL,
    strand varchar(10) NULL,
    assembly_id varchar(50) NULL,
    identity float8
); $$
AFTER LOAD DO
$$ DELETE FROM rnc_genome_mapping
WHERE identity = 1
AND assembly_id IN (SELECT distinct assembly_id from {table});
$$,
$$ INSERT INTO rnc_genome_mapping (rna_id, upi, taxid, region_id, chromosome,
                                   "start", stop, strand, assembly_id, "identity")
    (SELECT rna_id, upi, taxid, region_id, chromosome, "start", stop, strand,
            assembly_id, "identity"
	from {table})
; $$
;
"""


class GenomeMappingPGLoadExactMatches(PGLoader):  # pylint: disable=R0904
    """
    """
    assembly_id = luigi.Parameter()
    species = luigi.Parameter()
    taxid = luigi.IntParameter()
    division = luigi.Parameter()

    def control_filename(self):
        return os.path.join(genome_mapping().species(self.species),
                            'pgloader-exact-matches.ctl')

    def requires(self):
        return ParsePslOutput(taxid=self.taxid,
                              species=self.species,
                              assembly_id=self.assembly_id,
                              division=self.division)

    def control_file(self):
        filename = self.input()['exact'].path
        table = 'load_genome_mapping'
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table=table),
            table=table,
            search_path=self.db_search_path()
        )
