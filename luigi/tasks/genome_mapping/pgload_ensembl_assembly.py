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

import luigi

from tasks.utils.pgloader import PGLoader
from .update_ensembl_assembly import RetrieveEnsemblAssemblies
from .update_ensembl_assembly import RetrieveEnsemblGenomesAssemblies


CONTROL_FILE = """LOAD CSV
FROM '{filename}'
HAVING FIELDS
(
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    subdomain,
    example_chromosome,
    example_start,
    example_end,
    blat_mapping
)
INTO {db_url}
TARGET COLUMNS
(
    assembly_id,
    assembly_full_name,
    gca_accession,
    assembly_ucsc,
    common_name,
    taxid,
    ensembl_url,
    division,
    subdomain,
    example_chromosome,
    example_start,
    example_end,
    blat_mapping
)
WITH
    fields terminated by '\t'
SET
    search_path = 'rnacen'
BEFORE LOAD DO
$$ DROP TABLE IF EXISTS {table}; $$,
$$ CREATE TABLE rnacen.{table} (
	assembly_id varchar(255) NOT NULL,
	assembly_full_name varchar(255) NOT NULL,
	gca_accession varchar(20) NULL,
	assembly_ucsc varchar(100) NULL,
	common_name varchar(255) NOT NULL,
	taxid int4 NOT NULL,
	ensembl_url varchar(100) NULL,
	division varchar(20) NULL,
	blat_mapping int4 NULL,
	example_chromosome varchar(20) NULL,
	example_end int4 NULL,
	example_start int4 NULL,
	subdomain varchar(100) NOT NULL,
	CONSTRAINT {table}_pkey PRIMARY KEY (assembly_id)
); $$
AFTER LOAD DO
$$ DELETE FROM {permanent_table}
WHERE
assembly_id IN (SELECT assembly_id from {table});
$$,
$$ INSERT INTO {permanent_table} (
        assembly_id,
        assembly_full_name,
        gca_accession,
        assembly_ucsc,
        common_name,
        taxid,
        ensembl_url,
        division,
        blat_mapping,
        example_chromosome,
        example_end,
        example_start,
        subdomain
    )
    (SELECT
        assembly_id,
        assembly_full_name,
        gca_accession,
        assembly_ucsc,
        common_name,
        taxid,
        ensembl_url,
        division,
        blat_mapping,
        example_chromosome,
        example_end,
        example_start,
        subdomain
	FROM {table});
$$,
$$ DROP TABLE {table}; $$
;
"""


class GenomeMappingPGLoadEnsemblAssembly(PGLoader):  # pylint: disable=R0904
    """
    Update ensembl_assembly table with Ensembl data.
    """
    def requires(self):
        return RetrieveEnsemblAssemblies()

    def control_file(self):
        filename = RetrieveEnsemblAssemblies().output().fn
        table = 'load_ensembl_assembly'
        permanent_table = 'ensembl_assembly'
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table=table),
            table=table,
            permanent_table=permanent_table,
            search_path=self.db_search_path()
        )


class GenomeMappingPGLoadEnsemblGenomesAssembly(PGLoader):  # pylint: disable=R0904
    """
    Update ensembl_assembly table with Ensembl Genomes data.
    """
    def requires(self):
        return RetrieveEnsemblGenomesAssemblies()

    def control_file(self):
        filename = self.input().fn
        table = 'load_ensembl_assembly'
        permanent_table = 'ensembl_assembly'
        return CONTROL_FILE.format(
            filename=filename,
            db_url=self.db_url(table=table),
            table=table,
            permanent_table=permanent_table,
            search_path=self.db_search_path()
        )
